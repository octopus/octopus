!! Copyright (C) 2018 N. Tancogne-Dejean 
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

! ---------------------------------------------------------
!> This routine is an interface for constructing the orbital basis.
! ---------------------------------------------------------
subroutine X(orbitalbasis_build)(this, ions, mesh, kpt, ndim, skip_s_orb, use_all_orb, verbose)
  type(orbitalbasis_t), target, intent(inout)    :: this
  type(ions_t),         target, intent(in)       :: ions
  type(mesh_t),                 intent(in)       :: mesh
  type(distributed_t),          intent(in)       :: kpt
  integer,                      intent(in)       :: ndim
  logical,                      intent(in)       :: skip_s_orb 
  logical,                      intent(in)       :: use_all_orb
  logical, optional,            intent(in)       :: verbose

  integer :: ia, iorb, norb, offset, ios
  integer :: hubbardl, ii, nn, ll, mm, work, work2, iorbset
  FLOAT   :: hubbardj, jj
  integer :: n_s_orb
  type(orbitalset_t), pointer :: os
  logical :: hasjdependence
  logical :: verbose_, use_mesh

  PUSH_SUB(X(orbitalbasis_build))

  verbose_ = optional_default(verbose,.true.)

  !Do we use a mesh or a submesh to store the orbitals
  use_mesh = .not. this%submesh
#ifdef R_TCOMPLEX
  !In case of a phase, we want to use a submesh, to avoid problem
  !with periodic systems and self-overlapping submeshes
  use_mesh = .false.
#endif

  if(verbose_) then
    write(message(1),'(a)')    'Building the LDA+U localized orbital basis.'
    call messages_info(1)
  end if

  !We first count the number of orbital sets we have to treat
  norb = 0
  if( .not. use_all_orb ) then
    do ia = 1, ions%natoms
      hubbardl = species_hubbard_l(ions%atom(ia)%species)
      hubbardj = species_hubbard_j(ions%atom(ia)%species)
      if( hubbardl .eq. -1 ) cycle
  
      !This is a dirty way to detect if the pseudopotential has j-dependent atomic wavefunctions
      hasjdependence = .false.
      call species_iwf_j(ions%atom(ia)%species, 1, jj)
      if(abs(jj) > M_EPSILON) hasjdependence = .true.

      if(hasjdependence .and. abs(hubbardj) <= M_EPSILON) then
        norb = norb+2
      else
        norb = norb+1
      end if

      if( hubbardj /= 0 .and. .not. hasjdependence) then
        write(message(1),'(a,i1,a)') 'Atom ', ia, ' has no j-dependent atomic wavefunctions.'
        write(message(2),'(a)') 'This is not compatible with the hubbard_j option.'
        call messages_fatal(2)
      end if
    end do
  else
    do ia = 1, ions%natoms
      if(species_type(ions%atom(ia)%species) /= SPECIES_PSEUDO &
           .and. species_type(ions%atom(ia)%species) /= SPECIES_PSPIO) cycle
      work = 0
      n_s_orb = 0
      hubbardj = species_hubbard_j(ions%atom(ia)%species)
      do iorb = 1, species_niwfs(ions%atom(ia)%species)
        call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm) 
        call species_iwf_j(ions%atom(ia)%species, iorb, jj)
        if(ll == 0) n_s_orb = n_s_orb + 1
        work = max(work, ii)

        if( hubbardj /= 0 .and. abs(jj) > M_EPSILON ) then
          write(message(1),'(a,i1,a)') 'Atom ', ia, ' has no j-dependent atomic wavefunction.'
          write(message(2),'(a)') 'This is not compatible with the hubbard_j option.'
          call messages_fatal(2)  
        end if
      end do
      if(skip_s_orb) work = work-n_s_orb
      norb = norb + work
    end do
  end if

  if(norb == 0) then
    write(message(1),'(a)')  'No orbital set found. Please check your input file.'
    call messages_fatal(1)
  end if


  if(verbose_) then
    write(message(1),'(a, i3, a)')    'Found ', norb, ' orbital sets.'
    call messages_info(1)
  end if

  this%norbsets = norb
  SAFE_ALLOCATE(this%orbsets(1:norb))
  do iorbset = 1, this%norbsets
    call orbitalset_init(this%orbsets(iorbset))
  end do

  iorbset = 0
  do ia = 1, ions%natoms
    if(species_type(ions%atom(ia)%species) /= SPECIES_PSEUDO &
           .and. species_type(ions%atom(ia)%species) /= SPECIES_PSPIO) cycle

    hubbardj = species_hubbard_j(ions%atom(ia)%species)
    !This is a dirty way to detect if the pseudopotential has j-dependent atomic wavefunctions
    hasjdependence = .false.
    call species_iwf_j(ions%atom(ia)%species, 1, jj)
    if(abs(jj) >  M_EPSILON) hasjdependence = .true.
    if (debug%info .and. hasjdependence .and. verbose_) then
      write(message(1),'(a,i3,a)')  'Debug: Atom ', ia, ' has j-dependent pseudo-wavefunctions.'
      call messages_info(1)
    end if 

    if(.not. use_all_orb) then
      hubbardl = species_hubbard_l(ions%atom(ia)%species)
      if( hubbardl .eq. -1 ) cycle
      !In this case we only have one orbital or we only want one
      if(.not. hasjdependence .or. hubbardj /= M_ZERO &
            .or. hubbardl == 0 ) then
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(ions%atom(ia)%species)
          call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(ions%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(ions%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. hubbardj == jj ) then
            norb = norb + 1
            call orbitalset_set_jln(os, jj, hubbardl, nn)
            os%ii = ii
            os%radius = atomic_orbital_get_radius(ions, mesh, ia, iorb, 1, &
                          this%truncation, this%threshold)
          end if
        end do
        if( hasjdependence ) then
          os%ndim = ndim
          os%norbs = norb + int((os%jj - os%ll)*2)
        else
          os%ndim = 1
          os%norbs = norb
        end if
        os%Ueff = species_hubbard_u(ions%atom(ia)%species)
        os%alpha = species_hubbard_alpha(ions%atom(ia)%species)
        os%submesh = this%submesh
        os%spec => ions%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, ions, mesh, use_mesh, this%normalize)
      else
        !j = l-1/2
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(ions%atom(ia)%species)
          call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(ions%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(ions%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. jj < ll ) then
            norb = norb + 1
            call orbitalset_set_jln(os, jj, hubbardl, nn)
            os%ii = ii
            os%radius = atomic_orbital_get_radius(ions, mesh, ia, iorb, 1, &
                        this%truncation, this%threshold)
          end if
        end do
        os%ndim = ndim
        os%norbs = norb-1
        os%Ueff = species_hubbard_u(ions%atom(ia)%species)
        os%alpha = species_hubbard_alpha(ions%atom(ia)%species)
        os%submesh = this%submesh
        os%spec => ions%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, ions, mesh, use_mesh, this%normalize)

        !j = l+1/2
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(ions%atom(ia)%species)
          call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(ions%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(ions%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. jj > ll ) then
            norb = norb + 1
            call orbitalset_set_jln(os, jj, hubbardl, nn)
            os%ii = ii
            os%radius = atomic_orbital_get_radius(ions, mesh, ia, iorb, 1, &
                        this%truncation, this%threshold)
          end if
        end do
        os%ndim = ndim
        os%norbs = norb+1
        os%Ueff = species_hubbard_u(ions%atom(ia)%species)
        os%alpha = species_hubbard_alpha(ions%atom(ia)%species)
        os%submesh = this%submesh
        os%spec => ions%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, ions, mesh, use_mesh, this%normalize)
      end if
    else !use_all_orbitals
      ASSERT(.not.hasjdependence)
      work = 0
      n_s_orb = 0
      do iorb = 1, species_niwfs(ions%atom(ia)%species)
        call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
        if(ll == 0) n_s_orb = n_s_orb + 1
        work = max(work, ii)          
      end do
      offset = 0
      if(skip_s_orb) then
        work = work-n_s_orb
        offset = 1
      end if

      !We loop over the orbital sets of the atom ia
      do norb = 1, work
        os => this%orbsets(iorbset+norb)
        !We count the number of orbital for this orbital set
        work2 = 0
        do iorb = 1, species_niwfs(ions%atom(ia)%species)
          call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(ions%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(ions%atom(ia)%species, iorb, jj )
          if(ii == norb+offset .and. hubbardj == jj) then
            work2 = work2 + 1
            call orbitalset_set_jln(os, jj, ll, nn)
            os%ii = ii
            os%radius = atomic_orbital_get_radius(ions, mesh, ia, iorb, 1, &
                               this%truncation, this%threshold)
          end if
        end do
        os%norbs = work2
        os%ndim = 1
        os%Ueff = species_hubbard_u(ions%atom(ia)%species)
        os%alpha = species_hubbard_alpha(ions%atom(ia)%species)
        os%submesh = this%submesh
        os%spec => ions%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, ions, mesh, use_mesh, this%normalize)
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
    SAFE_ALLOCATE(os%phase(1:os%sphere%np, kpt%start:kpt%end))
    os%phase(:,:) = M_ZERO
    if(.not. this%submesh) then 
      SAFE_ALLOCATE(os%eorb_mesh(1:mesh%np, 1:os%norbs, 1:os%ndim, kpt%start:kpt%end))
      os%eorb_mesh(:,:,:,:) = M_ZERO
    else
      SAFE_ALLOCATE(os%eorb_submesh(1:os%sphere%np, 1:os%ndim, 1:os%norbs, kpt%start:kpt%end))
      os%eorb_submesh(:,:,:,:) = M_ZERO
    end if
  #endif

    ! We need to know the maximum number of points in order to allocate a temporary array
    ! to apply the phase in lda_u_apply
    if(os%sphere%np > this%max_np) this%max_np = os%sphere%np
  end do  

  do ios = 1, this%norbsets
    if(this%orbsets(ios)%sphere%np == -1) then
       write(message(1),'(a,a4,i1,a1,a)')    'Internal error: the orbital ',trim(species_label(this%orbsets(ios)%spec)), &
                      this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), ' has no grid point.'
       write(message(2),'(a)') 'Change the input file or use a pseudopotential that contains these orbitals.'
       call messages_fatal(2)
    end if
    if(verbose_) then
      write(message(1),'(a,i2,a,f8.5,a)')    'Orbital set ', ios, ' has a value of U of ',&
                         this%orbsets(ios)%Ueff   , ' Ha.'
      write(message(2),'(a,i2,a)')    'It contains ', this%orbsets(ios)%norbs, ' orbitals.'
      write(message(3),'(a,f8.5,a,i6,a)') 'The radius is ', this%orbsets(ios)%sphere%radius, &
                        ' Bohr,  with ', this%orbsets(ios)%sphere%np, ' grid points.'
       call messages_info(3)
    end if
  end do 

  this%size = 0
  do ios = 1, this%norbsets
    this%size = this%size + this%orbsets(ios)%norbs
  end do

  SAFE_ALLOCATE(this%global2os(2,this%size))
  SAFE_ALLOCATE(this%os2global(this%norbsets, this%maxnorbs))
  offset = 1
  do ios = 1, this%norbsets
    do iorb = 1, this%orbsets(ios)%norbs
      this%global2os(1,offset) = ios
      this%global2os(2,offset) = iorb
      this%os2global(ios, iorb) = offset
      offset = offset + 1
    end do
  end do

  POP_SUB(X(orbitalbasis_build))

end subroutine X(orbitalbasis_build)


! ---------------------------------------------------------
!> This routine constructd an empty orbital basis.
! ---------------------------------------------------------
subroutine X(orbitalbasis_build_empty)(this, mesh, kpt, ndim, nstates, verbose)
  type(orbitalbasis_t), target, intent(inout)    :: this
  type(distributed_t),          intent(in)       :: kpt
  type(mesh_t),         target, intent(in)       :: mesh
  integer,                      intent(in)       :: ndim
  integer,                      intent(in)       :: nstates
  logical, optional,            intent(in)       :: verbose

  integer :: ios, iorb, offset
  type(orbitalset_t), pointer :: os
  logical :: verbose_

  PUSH_SUB(X(orbitalbasis_build_empty))

  verbose_ = optional_default(verbose,.true.)

  if(verbose_) then
    write(message(1),'(a)')    'Building an empty LDA+U orbital basis.'
    call messages_info(1)
  end if

  ASSERT(nstates > 0)

  this%norbsets = 1
  SAFE_ALLOCATE(this%orbsets(1:1))
  call orbitalset_init(this%orbsets(1))
  os => this%orbsets(1)
  os%ii = -1
  os%radius = M_ZERO
  os%ndim = ndim
  os%norbs = nstates
  os%Ueff = M_ZERO
  os%alpha = M_ZERO
  os%submesh = .false.
  os%sphere%mesh => mesh
  nullify(os%spec)
  os%iatom = -1
  SAFE_ALLOCATE(os%X(orb)(1:mesh%np, 1:os%ndim, 1:os%norbs))
  os%X(orb)(:,:,:) = R_TOTYPE(M_ZERO)

  this%maxnorbs = nstates
  this%max_np = mesh%np 

  ! In case of complex wavefunction, we allocate the array for the phase correction
#ifdef R_TCOMPLEX
  SAFE_ALLOCATE(os%phase(1:mesh%np, kpt%start:kpt%end))
  os%phase(:,:) = M_ZERO
  SAFE_ALLOCATE(os%eorb_mesh(1:mesh%np, 1:os%norbs, 1:os%ndim, kpt%start:kpt%end))
  os%eorb_mesh(:,:,:,:) = M_ZERO
#endif

  this%size = nstates
  do ios = 1, this%norbsets
    this%size = this%size + this%orbsets(ios)%norbs
  end do

  SAFE_ALLOCATE(this%global2os(2,this%size))
  SAFE_ALLOCATE(this%os2global(this%norbsets, this%maxnorbs))
  offset = 1
  do ios = 1, this%norbsets
    do iorb = 1, this%orbsets(ios)%norbs
      this%global2os(1,offset) = ios
      this%global2os(2,offset) = iorb
      this%os2global(ios, iorb) = offset
      offset = offset + 1
    end do
  end do

  POP_SUB(X(orbitalbasis_build_empty))

end subroutine X(orbitalbasis_build_empty)
