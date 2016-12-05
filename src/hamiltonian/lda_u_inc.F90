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


subroutine X(hubbard_apply)(this, mesh, d, ik, psib, hpsib)
  implicit none 

  type(lda_u_t), intent(in)      :: this
  type(mesh_t),    intent(in)    :: mesh
  integer,         intent(in)    :: ik
  type(batch_t),   intent(in)    :: psib
  type(batch_t),   intent(inout) :: hpsib
  type(states_dim_t), intent(in) :: d

  integer :: ibatch, idim, ia, imp, im, ispin
  integer :: bind, is
  R_TYPE  :: weight, reduced
  R_TYPE, allocatable :: psi(:,:), hpsi(:,:)
  R_TYPE, allocatable :: dot(:)

  ! no push_sub because it is called too frequently

  SAFE_ALLOCATE(psi(1:mesh%np, 1:d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:d%dim))
  SAFE_ALLOCATE(dot(1:this%maxnorbs))

  ispin = states_dim_get_spin_index(d, ik)

  do ibatch = 1, psib%nst
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib,ibatch, mesh%np, hpsi)
    do ia = 1, this%natoms
      do idim = 1, d%dim
        ! We have to compute 
        ! hpsi> += sum_m |phi m> sum_m' Vmm' <phi m' | psi >
        !
        ! We first compute <phi m | psi> for all orbitals of the atom
        !
        do im = 1, this%norbs(ia)
          dot(im) = X(mf_dotp)(mesh, this%orbitals(im,ispin,ia)%X(orbital_mesh),& 
                                  psi(1:mesh%np,idim))
        end do
        !
        do im = 1, this%norbs(ia)
          ! sum_m' Vmm' <phi m' | psi >
          reduced = M_ZERO
          do imp = 1, this%norbs(ia)
            reduced = reduced + this%X(V)(im,imp,ispin,ia)*dot(imp)
          end do
        !
        !  call submesh_add_to_mesh(this%orbitals(im,ispin,ia)%sphere, &
          !                           this%orbitals(im,ispin,ia)%X(orbital), &
          !                           hpsi(:, idim), this%X(V)(im,imp,ispin,ia)*dot) 
  !        bind = batch_ist_idim_to_linear(psib, (/ibatch, idim/))
          !forall (is = 1:this%orbitals(im,ispin,ia)%sphere%np)
          !  hpsib%states_linear(bind)%X(psi)(this%orbitals(im,ispin,ia)%sphere%map(is)) = &
          !    hpsib%states_linear(bind)%X(psi)(this%orbitals(im,ispin,ia)%sphere%map(is)) &
          !         + this%orbitals(im,ispin,ia)%X(orbital)(is)*reduced
          !end forall

  !        forall (is = 1:mesh%np)
  !          hpsib%states_linear(bind)%X(psi)(is) = &
  !            hpsib%states_linear(bind)%X(psi)(is) &
  !                 + this%orbitals(im,ispin,ia)%X(orbital_mesh)(is)*reduced
  !        end forall
        !call submesh_add_to_mesh(this%orbitals(im,ispin,ia)%sphere, &
        !                         this%orbitals(im,ispin,ia)%X(orbital), &
        !                         hpsi(:, idim), reduced)
         call lalg_axpy(mesh%np, reduced, &
               this%orbitals(im,ispin,ia)%X(orbital_mesh), hpsi(1:mesh%np, idim))
       end do !im
     end do !idim
   end do !ia
   call batch_set_state(hpsib, ibatch, mesh%np, hpsi)
  end do !ibatch

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(dot)

end subroutine X(hubbard_apply)

! ---------------------------------------------------------
!> This routine compute the values of the occupation matrices
! ---------------------------------------------------------
subroutine X(update_occ_matrices)(this, geo, mesh, st, hubbard_dc)
  implicit none
 
  type(lda_u_t), intent(inout)    :: this
  type(geometry_t), intent(in)    :: geo
  type(mesh_t),     intent(in)    :: mesh
  type(states_t),  intent(inout)  :: st
  FLOAT, intent(inout)            :: hubbard_dc

  integer :: ia, im, ik, ist, ispin, norbs
  R_TYPE, allocatable :: psi(:,:)
  R_TYPE, allocatable :: dot(:)
  FLOAT   :: weight

  PUSH_SUB(update_occ_matrices)

  this%X(n)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:geo%natoms) = R_TOTYPE(M_ZERO)
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(dot(1:this%maxnorbs))

  !TODO: use symmetries of the occupation matrices
  do ik = st%d%kpt%start, st%d%kpt%end
    ispin =  states_dim_get_spin_index(st%d,ik)
    do ist = st%st_start, st%st_end
      weight = st%d%kweights(ik) * st%occ(ist, ik)

      call states_get_state(st, mesh, ist, ik, psi )   
      do ia = 1, geo%natoms
        norbs = this%norbs(ia)
        do im = 1, norbs
        !  dot(im) =  submesh_to_mesh_dotp(this%orbitals(im,ispin,ia)%sphere, st%d%dim, &
        !                        this%orbitals(im,ispin,ia)%X(orbital), psi)
          dot(im) = X(mf_dotp)(mesh, psi(1:mesh%np,1), this%orbitals(im,ispin,ia)%X(orbital_mesh))
        end do
        do im = 1, norbs
          this%X(n)(1:norbs,im,ispin,ia) = this%X(n)(1:norbs,im,ispin,ia) &
                                           + weight*dot(1:norbs)*R_CONJ(dot(im))
        end do
       !  call lalg_her( norbs, weight, this%X(n)(1:norbs,1:norbs,ispin,ia), dot)
      end do
    end do
  end do

  
  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(psi)

#if defined(HAVE_MPI)        
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp%comm, this%X(n))
  end if
#endif      

  call X(correct_energy_dc)(this, geo, st, hubbard_dc)
  call X(update_potential_lda_u)(this, geo, st)

  POP_SUB(update_occ_matrices)
end subroutine X(update_occ_matrices)

! ---------------------------------------------------------
!> This routine compute the value of the double counting term in the LDA+U energy
! ---------------------------------------------------------
subroutine X(correct_energy_dc)(this, geo, st, hubbard_dc)
  implicit none

  type(lda_u_t), intent(inout)    :: this
  type(geometry_t), intent(in)    :: geo
  type(states_t),  intent(in)     :: st
  FLOAT, intent(inout)         :: hubbard_dc

  integer :: ia, imp, im, ispin
  FLOAT :: U_I

  PUSH_SUB(correct_energy_dc)

  hubbard_dc = M_ZERO

  do ia = 1, geo%natoms
    U_I = species_hubbard_u(geo%atom(ia)%species)
    do ispin = 1, st%d%nspin
      !TODO: These are matrix operations, that could be optimized
      do im = 1, this%norbs(ia)
        do imp = 1, this%norbs(ia)
          hubbard_dc = hubbard_dc - CNST(0.5)*U_I*abs(this%X(n)(im,imp,ispin,ia))**2
        end do
        hubbard_dc = hubbard_dc + CNST(0.5)*U_I*this%X(n)(im,im,ispin,ia)
      end do
    end do
  end do

  POP_SUB(correct_energy_dc)
end subroutine X(correct_energy_dc)


! ---------------------------------------------------------
!> This routine compute the potential that, once multiplied
!> by the projector Pmm' and summed over m and m' for all the atoms
!> gives the full Hubbard potential
! ---------------------------------------------------------
subroutine X(update_potential_lda_u)(this, geo, st)
  implicit none

  type(lda_u_t), intent(inout)    :: this
  type(geometry_t), intent(in)    :: geo
  type(states_t),  intent(in)     :: st

  integer :: ia, im, ispin, norbs
  FLOAT :: U_I

  PUSH_SUB(update_potential_lda_u)

  !this%X(V)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:geo%natoms) = R_TOTYPE(M_ZERO)

  do ia = 1, geo%natoms
    U_I = species_hubbard_u(geo%atom(ia)%species)
    norbs = this%norbs(ia)
    do ispin = 1, st%d%nspin
      do im = 1, norbs
        this%X(V)(1:norbs,im,ispin,ia) = - U_I*this%X(n)(1:norbs,im,ispin,ia)
        this%X(V)(im,im,ispin,ia) = this%X(V)(im,im,ispin,ia) + CNST(0.5)*U_I
      end do
    end do
  end do

  POP_SUB(update_potential_lda_u)
end subroutine X(update_potential_lda_u)

! ---------------------------------------------------------
!> This routine computes [r,V_lda+u].
! ---------------------------------------------------------
!subroutine X(lda_u_commute_r)(this, mesh, ik, psi, gpsi)
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
!   PUSH_SUB(lda_u_commute_r)
!
!   SAFE_ALLOCATE(gspsi(1:mesh%np, 1:this%gs_st%d%dim))
!   SAFE_ALLOCATE(psi_r(1:mesh%np))
!   SAFE_ALLOCATE(tmpstate(1:mesh%np, 1:this%gs_st%d%dim))
!   
!   tmpstate(1:mesh%np, 1:this%gs_st%d%dim) = R_TOTYPE(M_ZERO)
!   do ia = 1, geo%natoms
!     do imp = 1, this%maxorb(ia)
!       dot = X(mf_dotp)(mesh, this%X(basis)(1:mesh%np,imp,ispin,ia),psi(1:mesh%np,idim))
!       do im = 1, this%maxorb(ia)
!          call lalg_axpy(mesh%np, -this%X(V)(im,imp,ispin,ia)*dot, &
!                 this%X(basis)(1:mesh%np,im,ispin,ia), hpsi(1:mesh%np, idim))
!         
!   do ist = 1, this%gs_st%nsti
!     call states_get_state(this%gs_st, mesh, ist, ik, gspsi )
!     !<gpsi|psi>
!     dot = X(mf_dotp)(mesh, this%gs_st%d%dim, gspsi(1:mesh%np,1:this%gs_st%d%dim), psi) &
!           * this%gs_st%occ(ist, ik)/ this%gs_st%smear%el_per_state
!     do idim = 1, this%gs_st%d%dim
!      call lalg_axpy(mesh%np, dot,  gspsi(1:mesh%np, idim), tmpstate(1:mesh%np, idim)) 
!     enddo
!   enddo
!   ! |gpsi> -= x|phim><phim'|psi>
!   do idim = 1, this%gs_st%d%dim
!     do idir = 1, mesh%sb%dim
!       gpsi(1:mesh%np, idir, idim) = gpsi(1:mesh%np, idir, idim) &
!             -  mesh%x(1:mesh%np, idir) * tmpstate(1:mesh%np, idim)
!     enddo
!   enddo
!
!   do idim = 1, this%gs_st%d%dim
!     do idir = 1, mesh%sb%dim
!  
!       psi_r(1:mesh%np) = mesh%x(1:mesh%np,idir) * psi(1:mesh%np, idim)
!
!       tmpstate(1:mesh%np,idim) = R_TOTYPE(M_ZERO)
!       do ist = 1, this%gs_st%nst
!         call states_get_state(this%gs_st, mesh, ist, ik, gspsi )
!         ! <gspsi|r|psi>
!         dot = X(mf_dotp)(mesh, gspsi(1:mesh%np,idim), psi_r(1:mesh%np)) &
!           * this%gs_st%occ(ist, ik) / this%gs_st%smear%el_per_state
!         call lalg_axpy(mesh%np, dot,  gspsi(1:mesh%np, idim), tmpstate(1:mesh%np, idim))
!       enddo
!       ! |gpsi> +=  |gspsi><gspsi|x|psi>
!       gpsi(1:mesh%np, idir, idim) = gpsi(1:mesh%np, idir, idim) &
!             +  tmpstate(1:mesh%np, idim)
!     enddo
!   enddo
!
!   SAFE_DEALLOCATE_A(gspsi)
!   SAFE_DEALLOCATE_A(psi_r)
!   SAFE_DEALLOCATE_A(tmpstate) 
!
!   POP_SUB(lda_u_commute_r)
!end subroutine X(lda_u_commute_r)

! ---------------------------------------------------------
!> This routine is an interface for constructing the orbital basis.
! ---------------------------------------------------------
subroutine X(construct_orbital_basis)(this, geo, mesh, st)

  implicit none

  type(lda_u_t),             intent(inout)    :: this
  type(geometry_t), target,  intent(in)       :: geo
  type(mesh_t),              intent(in)       :: mesh
  type(states_t),            intent(in)       :: st 

  integer :: ia, iorb, norb, ispin
  integer ::  hubbardl, ii, ll, mm
  FLOAT   :: norm !, norm_proj


  PUSH_SUB(X(construct_orbital_basis))

  write(message(1),'(a)')    'Building the LDA+U localized orbital basis.'
  call messages_info(1)

  !We first find the number of orbitals for each atoms
  SAFE_ALLOCATE(this%norbs(1:geo%natoms))
  this%norbs = 0
  do ia = 1, geo%natoms
    hubbardl = species_hubbard_l(geo%atom(ia)%species)
    if( hubbardl .eq. 0 ) cycle
    do ispin = 1,st%d%nspin
      this%norbs(ia) = 0
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, ispin, ii, ll, mm)
        if(ll .eq. hubbardl ) this%norbs(ia) = this%norbs(ia) + 1
      end do !iorb
    end do !ispin
  end do  !ia
  this%maxnorbs = maxval(this%norbs) 
  
  do ia = 1, geo%natoms
    write(message(1),'(a,i2,a,f8.5,a)')    'Atom ', ia, ' has a value of U of ',&
                            species_hubbard_u(geo%atom(ia)%species), ' Ha.'
    call messages_info(1)

    write(message(1),'(a,i2,a, i3)')    'Found ', this%norbs(ia), ' orbitals for atom ', ia
     call messages_info(1)
  end do
 
  SAFE_ALLOCATE(this%orbitals(1:this%maxnorbs, 1:st%d%nspin, 1:geo%natoms))

  do ia = 1, geo%natoms
    hubbardl = species_hubbard_l(geo%atom(ia)%species)
    if( hubbardl .eq. M_ZERO ) cycle
 
    do ispin = 1,st%d%nspin
      norb = 0
!      norm_proj = CNST(0.0)
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, ispin, ii, ll, mm)
        if(ll .eq. hubbardl ) then 
          norb = norb + 1
          !We obtain the orbital
          call X(get_atomic_orbital)(geo, mesh, ia, iorb, ispin, this%orbitals(norb,ispin,ia),&
                                     this%truncate) 
          ! We have to normalize the orbitals, 
          ! in case the orbitals that comes out of the pseudo are not properly normalised
          norm = X(mf_nrm2)(mesh, this%orbitals(norb,ispin,ia)%X(orbital_mesh)(1:mesh%np))
          this%orbitals(norb,ispin,ia)%X(orbital_mesh)(1:mesh%np) =  &
                  this%orbitals(norb,ispin,ia)%X(orbital_mesh)(1:mesh%np) /sqrt(norm)
         !   norm_proj = norm_proj + norm
        endif
      end do
      
!      !Normalizing the projector
!      do iorb = 1, norb
!        this%X(basis)(1:mesh%np, ispin, iorb, ia) = this%X(basis)(1:mesh%np, ispin, iorb, ia) / sqrt(norm_proj)
!      end do
    end do
  end do

  POP_SUB(X(construct_orbital_basis))

end subroutine X(construct_orbital_basis)


! ---------------------------------------------------------
!> This routine returns the atomic orbital basis -- provided
!! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(get_atomic_orbital) (geo, mesh, iatom, iorb, ispin, orb, truncate)
  type(mesh_t),             intent(in)    :: mesh
  type(geometry_t), target, intent(in)    :: geo
  integer,                  intent(in)    :: iatom
  integer,                  intent(in)    :: iorb
  integer,                  intent(in)    :: ispin
  type(orbital_t),          intent(inout) :: orb
  logical,                  intent(in)    :: truncate

  type(species_t), pointer :: spec
  integer :: ii, ll, mm, ispin_
  FLOAT :: radius
  logical :: complex_ylms

  PUSH_SUB(X(get_atomic_orbital))

  ispin_ = ispin 

  spec => geo%atom(iatom)%species
  ASSERT(iorb <= species_niwfs(spec))

  nullify(orb%dorbital_sphere)
  nullify(orb%zorbital_sphere)
  nullify(orb%dorbital_mesh)
  nullify(orb%zorbital_mesh)

  call species_iwf_ilm(spec, iorb, ispin_, ii, ll, mm)

  radius = species_get_iwf_radius(spec, ii, ispin_) 

  !If asked, we truncate the orbital to the radius on the projector spheres 
  !of the NL part of the pseudopotential.
  !This is a way to garanty no overlap between orbitals of different atoms.
  if(truncate .and. species_is_ps(spec)) &
    radius = min(radius,species_get_ps_radius(spec))

  ! make sure that if the spacing is too large, the orbitals fit in a few points at least
  radius = max(radius, CNST(2.0)*maxval(mesh%spacing(1:mesh%sb%dim)))

 
  !We initialise the submesh corresponding to the orbital 
  call submesh_init(orb%sphere, mesh%sb, mesh, geo%atom(iatom)%x, radius)

  !We allocate both the orbital on the submesh and on the complete mesh
  SAFE_ALLOCATE(orb%X(orbital_sphere)(1:orb%sphere%np))
  SAFE_ALLOCATE(orb%X(orbital_mesh)(1:mesh%np))
  orb%X(orbital_sphere) = M_ZERO
  orb%X(orbital_mesh) = M_ZERO

  !This is a bit dirty. This is the default behavior, see LCAO.
  complex_ylms = .false.

  !We get the orbital from the pseudopotential
  #ifdef R_TCOMPLEX
  if(.not. complex_ylms) then
    !In this case we want to get a real orbital and to store it in complex array
    SAFE_ALLOCATE(orb%dorbital_sphere(1:orb%sphere%np))
    call dspecies_get_orbital_submesh(spec, orb%sphere, ii, ll, mm, ispin_, geo%atom(iatom)%x, orb%dorbital_sphere)
    call submesh_add_to_mesh(orb%sphere, orb%dorbital_sphere, orb%X(orbital_mesh))
    orb%X(orbital_sphere)(1:orb%sphere%np) = orb%dorbital_sphere(1:orb%sphere%np)
    SAFE_DEALLOCATE_P(orb%dorbital_sphere)
  else
  #endif
    call X(species_get_orbital_submesh)(spec, orb%sphere, ii, ll, mm, ispin_, geo%atom(iatom)%x, orb%X(orbital_sphere))
    call submesh_add_to_mesh(orb%sphere, orb%X(orbital_sphere), orb%X(orbital_mesh))
  #ifdef R_TCOMPLEX
  end if
  #endif

  POP_SUB(X(get_atomic_orbital))

end subroutine X(get_atomic_orbital)

! ---------------------------------------------------------

