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
!> This routine fills state psi with an atomic orbital -- provided
!! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(lcao_atomic_orbital) (this, iorb, mesh, st, geo, psi, spin_channel)
  type(lcao_t),             intent(in)    :: this
  integer,                  intent(in)    :: iorb
  type(mesh_t),             intent(in)    :: mesh
  type(states_t),           intent(in)    :: st
  type(geometry_t), target, intent(in)    :: geo
  R_TYPE,                   intent(inout) :: psi(:, :)
  integer,                  intent(in)    :: spin_channel

  type(species_t), pointer :: spec
  integer :: idim, iatom, jj, ip
  FLOAT, allocatable :: ao(:)
  type(profile_t), save :: prof

  call profiling_in(prof, "ATOMIC_ORBITAL")
  PUSH_SUB(X(lcao_atomic_orbital))

  ASSERT(iorb >= 1)
  ASSERT(iorb <= this%maxorbs)

  psi(1:mesh%np, 1:st%d%dim) = R_TOTYPE(M_ZERO)

  iatom = this%atom(iorb)
  jj = this%level(iorb)
  idim = this%ddim(iorb)
  spec => geo%atom(iatom)%spec
  ASSERT(jj <= species_niwfs(spec))

  SAFE_ALLOCATE(ao(1:mesh%np))

  call species_get_orbital(spec, mesh, jj, max(spin_channel, idim), geo%atom(iatom)%x, ao)
  
  do ip = 1, mesh%np
    psi(ip, idim) = ao(ip)
  end do

  SAFE_DEALLOCATE_A(ao)

  POP_SUB(X(lcao_atomic_orbital))
  call profiling_out(prof)

end subroutine X(lcao_atomic_orbital)


! ---------------------------------------------------------
subroutine X(lcao_wf)(this, st, gr, geo, hm, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer, optional,   intent(in)    :: start

  integer :: nst, ik, n1, n2, idim, lcao_start, ie, max
  R_TYPE, allocatable :: hpsi(:, :, :), overlap(:, :, :)
  FLOAT, allocatable :: ev(:)
  R_TYPE, allocatable :: hamilt(:, :, :), lcaopsi(:, :, :), lcaopsi2(:, :)
  integer :: kstart, kend, ispin, iunit_h, iunit_s, iunit_e

#ifdef HAVE_MPI
  FLOAT, allocatable :: tmp(:, :)
#endif

  PUSH_SUB(X(lcao_wf))
  
  write(message(1),'(a,i6,a)') 'Info: Performing initial LCAO calculation with ', &
       this%norbs,' orbitals.'
  call messages_info(1)
          

  nst = st%nst
  kstart = st%d%kpt%start
  kend = st%d%kpt%end

  lcao_start = optional_default(start, 1)
  if(st%parallel_in_states .and. st%st_start > lcao_start) lcao_start = st%st_start

  ! Allocation of variables

  SAFE_ALLOCATE(lcaopsi(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%spin_channels))
  SAFE_ALLOCATE(lcaopsi2(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(hpsi(gr%mesh%np, st%d%dim, kstart:kend))
  SAFE_ALLOCATE(hamilt(this%norbs, this%norbs, kstart:kend))
  SAFE_ALLOCATE(overlap(this%norbs, this%norbs, st%d%spin_channels))

  ie = 0
  max = this%norbs * (this%norbs + 1)/2

  message(1) = "Info: Getting Hamiltonian matrix elements."
  call messages_info(1)

  if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max)

#ifdef LCAO_DEBUG
! This code (and related below) is commented out because it causes mysterious optimization
! problems with PGI 12.4.0 -- LAPACK fails in diagonalization, only if mesh partition from scratch!
  if(this%write_matrices) then
    iunit_h = io_open(file=trim(STATIC_DIR)//'lcao_hamiltonian', action='write')
    iunit_s = io_open(file=trim(STATIC_DIR)//'lcao_overlap', action='write')
    iunit_e = io_open(file=trim(STATIC_DIR)//'lcao_eigenvectors', action='write')
    if(mpi_grp_is_root(mpi_world)) then
      write(iunit_h,'(4a6,a15)') 'iorb', 'jorb', 'ik', 'spin', 'hamiltonian'
      write(iunit_s,'(3a6,a15)') 'iorb', 'jorb', 'spin', 'overlap'
      write(iunit_e,'(4a6,a15)') 'ieig', 'jorb', 'ik', 'spin', 'coefficient'
    endif
  endif
#endif

  do n1 = 1, this%norbs
    
    do ispin = 1, st%d%spin_channels
      call X(get_ao)(this, st, gr%mesh, geo, n1, ispin, lcaopsi(:, :, ispin), use_psi = .true.)
    end do

    do ik = kstart, kend
      ispin = states_dim_get_spin_index(st%d, ik)
      call X(hamiltonian_apply)(hm, gr%der, lcaopsi(:, :, ispin), hpsi(:, :, ik), n1, ik)
    end do

    do n2 = n1, this%norbs
      do ispin = 1, st%d%spin_channels

        call X(get_ao)(this, st, gr%mesh, geo, n2, ispin, lcaopsi2, use_psi = .true.)

        overlap(n1, n2, ispin) = X(mf_dotp)(gr%mesh, st%d%dim, lcaopsi(:, :, ispin), lcaopsi2)
        overlap(n2, n1, ispin) = R_CONJ(overlap(n1, n2, ispin))
#ifdef LCAO_DEBUG
        if(this%write_matrices .and. mpi_grp_is_root(mpi_world)) &
          write(iunit_s,'(3i6,2f15.6)') n1, n2, ispin, overlap(n1, n2, ispin)
#endif

        do ik = kstart, kend
          if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle
          hamilt(n1, n2, ik) = X(mf_dotp)(gr%mesh, st%d%dim, hpsi(:, :, ik), lcaopsi2)
          hamilt(n2, n1, ik) = R_CONJ(hamilt(n1, n2, ik))

#ifdef LCAO_DEBUG
          if(this%write_matrices .and. mpi_grp_is_root(mpi_world)) &
            write(iunit_h,'(4i6,2f15.6)') n1, n2, ik, ispin, units_from_atomic(units_out%energy, hamilt(n1, n2, ik))
#endif
        end do
      end do
      
      ie = ie + 1
    end do

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ie, max)
  end do

  if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

  SAFE_DEALLOCATE_A(hpsi)

  SAFE_ALLOCATE(ev(1:this%norbs))

  do ik = kstart, kend
    ispin = states_dim_get_spin_index(st%d, ik)
    call lalg_geneigensolve(this%norbs, hamilt(:, :, ik), overlap(:, :, ispin), ev)

#ifdef HAVE_MPI
    ! the eigenvectors are not unique due to phases and degenerate subspaces, but
    ! they must be consistent among processors in domain parallelization
    if(gr%mesh%parallel_in_domains) &
      call MPI_Bcast(hamilt(1, 1, ik), this%norbs**2, R_MPITYPE, 0, gr%mesh%mpi_grp%comm, mpi_err)
#endif

    st%eigenval(lcao_start:nst, ik) = ev(lcao_start:nst)

    if(lcao_start /= 1) then 
      ASSERT(associated(st%X(psi)))
      st%X(psi)(1:gr%mesh%np, 1:st%d%dim, lcao_start:st%st_end, ik) = R_TOTYPE(M_ZERO)
    end if
  end do

#ifdef LCAO_DEBUG
  if(this%write_matrices .and. mpi_grp_is_root(mpi_world)) then
    do ik =  kstart, kend
      do n2 = 1, this%norbs
        do n1 = 1, this%norbs
          write(iunit_e,'(4i6,2f15.6)') n2, n1, ik, ispin, hamilt(n1, n2, ik)
        enddo
      enddo
    enddo

    call io_close(iunit_h)
    call io_close(iunit_s)
    call io_close(iunit_e)
  endif
#endif

#ifdef HAVE_MPI
  if(st%d%kpt%parallel) then
    ASSERT(.not. st%parallel_in_states)
    SAFE_ALLOCATE(tmp(1:st%nst, kstart:kend))
    tmp(1:st%nst, kstart:kend) = st%eigenval(1:st%nst, kstart:kend)
    call MPI_Allgatherv(tmp(:, kstart:), st%nst * (kend - kstart + 1), MPI_FLOAT, &
         st%eigenval, st%d%kpt%num(:) * st%nst, (st%d%kpt%range(1, :) - 1) * st%nst, MPI_FLOAT, &
         st%d%kpt%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(tmp)
  end if
#endif

  if(lcao_start == 1) call states_set_zero(st)

  ! Change of basis
  do n2 = 1, this%norbs
    do ispin = 1, st%d%spin_channels
      
      call X(get_ao)(this, st, gr%mesh, geo, n2, ispin, lcaopsi2, use_psi = .false.)

      do ik =  kstart, kend
        if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle
        do idim = 1, st%d%dim
          do n1 = lcao_start, st%st_end
            call states_get_state(st, gr%mesh, idim, n1, ik, lcaopsi(:, 1, 1))
            call lalg_axpy(gr%mesh%np, hamilt(n2, n1, ik), lcaopsi2(:, idim), lcaopsi(:, 1, 1))
            call states_set_state(st, gr%mesh, idim, n1, ik, lcaopsi(:, 1, 1))
          end do
        end do
      end do
      
    end do
  end do

  SAFE_DEALLOCATE_A(ev)
  SAFE_DEALLOCATE_A(hamilt)
  SAFE_DEALLOCATE_A(overlap)
  
  SAFE_DEALLOCATE_A(lcaopsi)
  SAFE_DEALLOCATE_A(lcaopsi2)

  POP_SUB(X(lcao_wf))
end subroutine X(lcao_wf)

! ---------------------------------------------------------
subroutine X(init_orbitals)(this, st, gr, geo, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  integer, optional,   intent(in)    :: start

  integer :: iorb, ispin, ist, ik, size
  integer :: nst, kstart, kend, lcao_start
  R_TYPE, allocatable :: ao(:, :)

  PUSH_SUB(X(init_orbitals))

  nst = st%nst
  kstart = st%d%kpt%start
  kend = st%d%kpt%end

  lcao_start = optional_default(start, 1)
  if(st%parallel_in_states .and. st%st_start > lcao_start) lcao_start = st%st_start

  ! We calculate the atomic orbitals first. To save memory we put
  ! all the orbitals we can in the part of st%Xpsi that we are going
  ! to overwrite and then the rest is stored in a single-precision
  ! buffer.

  SAFE_ALLOCATE(this%cst(1:this%norbs, 1:st%d%spin_channels))
  SAFE_ALLOCATE(this%ck(1:this%norbs, 1:st%d%spin_channels))
  SAFE_ALLOCATE(ao(1:gr%mesh%np, 1:st%d%dim))

  this%ck = 0

  iorb = 1
  ispin = 1

  ! first store in st%Xpsi
  ist_loop: do ist = lcao_start, st%st_end
    do ik = kstart, kend

      this%cst(iorb, ispin) = ist
      this%ck(iorb, ispin) = ik

      call X(lcao_atomic_orbital)(this, iorb, gr%mesh, st, geo, ao, ispin)
      call states_set_state(st, gr%mesh, ist, ik, ao)

      if(ispin < st%d%spin_channels) then
        ispin = ispin + 1
      else
        ispin = 1
        iorb = iorb + 1
        if(iorb > this%norbs) exit ist_loop
      end if

    end do
  end do ist_loop

  if(ispin < st%d%spin_channels) iorb = iorb - 1 ! we have not completed all the spin channels

  ! if there are any orbitals left, allocate extra space for them

  if(iorb <= this%norbs) then

    size = (this%norbs - iorb + 1) * st%d%spin_channels
    write(message(1), '(a, i5, a)') "Info: Single-precision storage for ", size, " extra orbitals will be allocated."
    call messages_info(1)

    SAFE_ALLOCATE(this%X(buff)(1:gr%mesh%np, 1:st%d%dim, iorb:this%norbs, 1:st%d%spin_channels))

    do iorb = iorb, this%norbs
      do ispin = 1, st%d%spin_channels
        call X(lcao_atomic_orbital)(this, iorb, gr%mesh, st, geo, ao, ispin)
        this%X(buff)(1:gr%mesh%np, 1:st%d%dim, iorb, ispin) = ao(1:gr%mesh%np, 1:st%d%dim)
      end do
    end do

  end if

  SAFE_DEALLOCATE_A(ao)

  POP_SUB(X(init_orbitals))

end subroutine X(init_orbitals)


! ---------------------------------------------------------
subroutine X(get_ao)(this, st, mesh, geo, iorb, ispin, ao, use_psi)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(in)    :: st
  type(mesh_t),        intent(in)    :: mesh
  type(geometry_t),    intent(in)    :: geo
  integer,             intent(in)    :: iorb
  integer,             intent(in)    :: ispin
  R_TYPE,              intent(out)   :: ao(:, :)
  logical,             intent(in)    :: use_psi
  
  PUSH_SUB(X(get_ao))
  
  if(this%ck(iorb, ispin) == 0) then
    ao(1:mesh%np, 1:st%d%dim) = this%X(buff)(1:mesh%np, 1:st%d%dim, iorb, ispin)
  else
    if(use_psi) then
      call states_get_state(st, mesh, this%cst(iorb, ispin), this%ck(iorb, ispin), ao)
    else
      call X(lcao_atomic_orbital)(this, iorb, mesh, st, geo, ao, ispin)
    end if
  end if

  POP_SUB(X(get_ao))
  
end subroutine X(get_ao)

! ---------------------------------------------------------

subroutine X(lcao_alt_init_orbitals)(this, st, gr, geo, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  integer, optional,   intent(in)    :: start

  integer :: iatom, norbs, dof

  PUSH_SUB(X(lcao_alt_init_orbitals))

  if(present(start)) then
    ASSERT(start == 1)
  end if

  write(message(1), '(a,i6,a)') 'Info: Performing LCAO calculation with ', this%nbasis, ' orbitals.'
  write(message(2), '(a)') ' '
  call messages_info(2)

  if (this%nbasis < st%nst) then
    write(message(1), '(a)') 'Not enough atomic orbitals to initialize all states,'
    write(message(2), '(i6,a)') st%nst - this%nbasis, ' states will be randomized.'
    if(this%derivative) then
      call messages_warning(2)
    else
      write(message(3), '(a)') 'You can double the number of atomic orbitals by setting'
      write(message(4), '(a)') 'LCAOExtraOrbitals to yes.'
      call messages_warning(4)
    end if
  end if

  SAFE_ALLOCATE(this%sphere(1:geo%natoms))
  SAFE_ALLOCATE(this%orbitals(1:geo%natoms))

  dof = 0
  do iatom = 1, geo%natoms
    norbs = species_niwfs(geo%atom(iatom)%spec)

    ! initialize the radial grid
    call submesh_init_sphere(this%sphere(iatom), gr%mesh%sb, gr%mesh, geo%atom(iatom)%x, this%radius(iatom))
    INCR(dof, this%sphere(iatom)%np*this%mult*norbs)
    call batch_init(this%orbitals(iatom), 1, this%mult*norbs)
  end do

  if(this%keep_orb) then
    write(message(1), '(a,f10.2,a)') &
      'Info: LCAO requires', dof*CNST(8.0)/CNST(1024.0)**2, ' Mb of memory for atomic orbitals.'
    call messages_info(1)
  end if

  POP_SUB(X(lcao_alt_init_orbitals))
end subroutine X(lcao_alt_init_orbitals)

! ---------------------------------------------------------
!> The alternative implementation.
subroutine X(lcao_alt_wf) (this, st, gr, geo, hm, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: start

  integer :: iatom, jatom, ik, ist, ispin, nev, ib
  integer :: ibasis, jbasis, iorb, jorb, norbs
  R_TYPE, allocatable :: hamiltonian(:, :), overlap(:, :), aa(:, :), bb(:, :)
  integer :: prow, pcol, ilbasis, jlbasis
  R_TYPE, allocatable :: psii(:, :, :), hpsi(:, :, :)
  type(batch_t) :: hpsib, psib
  FLOAT, allocatable :: eval(:)
  R_TYPE, allocatable :: evec(:, :), levec(:, :)  
  FLOAT :: dist2
  type(profile_t), save :: prof_matrix, prof_wavefunction

  PUSH_SUB(X(lcao_alt_wf))

  ASSERT(start == 1)

  if(.not. this%parallel) then
    SAFE_ALLOCATE(hamiltonian(1:this%nbasis, 1:this%nbasis))
    SAFE_ALLOCATE(overlap(1:this%nbasis, 1:this%nbasis))
  else
    SAFE_ALLOCATE(hamiltonian(1:this%lsize(1), 1:this%lsize(2)))
    SAFE_ALLOCATE(overlap(1:this%lsize(1), 1:this%lsize(2)))
  end if

  SAFE_ALLOCATE(aa(1:this%maxorb, 1:this%maxorb))
  SAFE_ALLOCATE(bb(1:this%maxorb, 1:this%maxorb))  
  SAFE_ALLOCATE(psii(1:gr%mesh%np_part, 1:st%d%dim, this%maxorb))
  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim, this%maxorb))

  call states_set_zero(st)

  do ispin = 1, st%d%spin_channels

    if(st%d%spin_channels > 1) then
      write(message(1), '(a)') ' '
      write(message(2), '(a,i1)') 'LCAO for spin channel ', ispin
      write(message(3), '(a)') ' '
      call messages_info(3)
    end if

    if(ispin > 1) then
      ! we need to deallocate previous orbitals
      do iatom = 1, geo%natoms
        call lcao_alt_end_orbital(this%orbitals(iatom))
      end do
    end if

    ! iterate over the kpoints for this spin
    do ik = st%d%kpt%start, st%d%kpt%end
      if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle

      if(st%d%nik > st%d%spin_channels) then
        write(message(1), '(a)') ' '
        write(message(2), '(a,i5)') 'LCAO for k-point ', states_dim_get_kpoint_index(st%d, ik)
        write(message(3), '(a)') ' '
        call messages_info(3)
      end if

      message(1) = 'Calculating matrix elements.'
      call messages_info(1)

      call profiling_in(prof_matrix, "LCAO_MATRIX")

      hamiltonian = R_TOTYPE(M_ZERO)
      overlap = R_TOTYPE(M_ZERO)

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, geo%natoms)

      do iatom = 1, geo%natoms
        norbs = this%norb_atom(iatom)

        call lcao_alt_get_orbital(this%orbitals(iatom), this%sphere(iatom), geo, ispin, iatom, this%norb_atom(iatom))

        psii = M_ZERO

        call batch_init( psib, st%d%dim, this%atom_orb_basis(iatom, 1), this%atom_orb_basis(iatom, norbs), psii)
        call batch_init(hpsib, st%d%dim, this%atom_orb_basis(iatom, 1), this%atom_orb_basis(iatom, norbs), hpsi)

        call X(submesh_batch_add)(this%sphere(iatom), this%orbitals(iatom), psib)
        call X(hamiltonian_apply_batch)(hm, gr%der, psib, hpsib, ik)

        do jatom = 1, geo%natoms
          if(.not. this%calc_atom(jatom)) cycle
          ! we only calculate the upper triangle
          if(jatom < iatom) cycle

          dist2 = sum((geo%atom(iatom)%x(1:MAX_DIM) - geo%atom(jatom)%x(1:MAX_DIM))**2)

          if(dist2 > (this%radius(iatom) + this%radius(jatom) + this%lapdist)**2) cycle

          call lcao_alt_get_orbital(this%orbitals(jatom), this%sphere(jatom), geo, ispin, jatom, this%norb_atom(jatom))

          ibasis = this%atom_orb_basis(iatom, 1)
          jbasis = this%atom_orb_basis(jatom, 1)

          call X(submesh_batch_dotp_matrix)(this%sphere(jatom), hpsib, this%orbitals(jatom), aa)

          if(dist2 > (this%radius(iatom) + this%radius(jatom))**2) then 
            bb = M_ZERO
          else
            call X(submesh_batch_dotp_matrix)(this%sphere(jatom), psib, this%orbitals(jatom), bb)
          end if

          if(.not. this%keep_orb) call lcao_alt_end_orbital(this%orbitals(jatom))

          !now, store the result in the matrix

          if(.not. this%parallel) then
            do iorb = 1, norbs
              do jorb = 1, this%norb_atom(jatom)
                hamiltonian(ibasis - 1 + iorb, jbasis - 1 + jorb) = aa(iorb, jorb)
                overlap(ibasis - 1 + iorb, jbasis - 1 + jorb) = bb(iorb, jorb)
              end do
            end do
          else
            do iorb = 1, norbs
              do jorb = 1, this%norb_atom(jatom)

                call lcao_local_index(this, ibasis - 1 + iorb,  jbasis - 1 + jorb, &
                  ilbasis, jlbasis, prow, pcol)
                if(all((/prow, pcol/) == this%myroc)) then
                  hamiltonian(ilbasis, jlbasis) = aa(iorb, jorb)
                  overlap(ilbasis, jlbasis) = bb(iorb, jorb)
                end if
              end do
            end do
          end if

        end do ! jatom

        call batch_end(psib)
        call batch_end(hpsib)

        if(.not. this%keep_orb) call lcao_alt_end_orbital(this%orbitals(iatom))

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(iatom, geo%natoms)
      end do

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

      message(1) = 'Diagonalizing Hamiltonian.'
      call messages_info(1)

      call profiling_out(prof_matrix)

      ! the number of eigenvectors we need
      nev = min(this%nbasis, st%nst) 

      SAFE_ALLOCATE(eval(1:this%nbasis))

      if(this%parallel) then
        SAFE_ALLOCATE(levec(1:this%lsize(1), 1:this%lsize(2)))
        SAFE_ALLOCATE(evec(1:this%nbasis, st%st_start:st%st_end))
      else 
        SAFE_ALLOCATE(evec(1:this%nbasis, 1:this%nbasis))
      end if

      call diagonalization()

      call profiling_in(prof_wavefunction, "LCAO_WAVEFUNCTIONS")

      message(1) = 'Generating wavefunctions.'
      call messages_info(1)

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, this%nbasis)

      ! set the eigenvalues
      st%eigenval(1:nev, ik) = eval(1:nev)
      st%eigenval(this%nbasis + 1:st%nst, ik) = HUGE(M_ONE)

      ibasis = 1
      do iatom = 1, geo%natoms
        norbs = this%norb_atom(iatom)

        call lcao_alt_get_orbital(this%orbitals(iatom), this%sphere(iatom), geo, ispin, iatom, this%norb_atom(iatom))

        do ib = st%block_start, st%block_end
          call X(submesh_batch_add_matrix)(this%sphere(iatom), evec(ibasis:, states_block_min(st, ib):), &
            this%orbitals(iatom), st%psib(ib, ik))
        end do

        if(.not. this%keep_orb) call lcao_alt_end_orbital(this%orbitals(iatom))

        ibasis = ibasis + norbs
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ibasis, this%nbasis)
      end do

      SAFE_DEALLOCATE_A(eval)
      SAFE_DEALLOCATE_A(evec)
      SAFE_DEALLOCATE_A(levec)

      ! if we do not have enough basis functions randomize the missing orbitals
      do ist = this%nbasis + 1, st%st_end
        call X(mf_random)(gr%mesh, psii(:, 1, 1))
        call states_set_state(st, gr%mesh, 1, ist, ik, psii(:, 1, 1))
      end do

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')
      call profiling_out(prof_wavefunction)
    end do
  end do

  do iatom = 1, geo%natoms
    call submesh_end(this%sphere(iatom))
    call lcao_alt_end_orbital(this%orbitals(iatom))
    call batch_end(this%orbitals(iatom))
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(hpsi)  
  SAFE_DEALLOCATE_A(hamiltonian)
  SAFE_DEALLOCATE_A(overlap)

  POP_SUB(X(lcao_alt_wf))

contains

  subroutine diagonalization()
    integer              :: neval_found, info, lwork
    R_TYPE               :: tmp(3) !< size must be at least 3 according to ScaLAPACK
    R_TYPE,  allocatable :: work(:)
    integer, allocatable :: iwork(:), ifail(:)
#ifdef R_TCOMPLEX
    FLOAT,   allocatable :: rwork(:)
#endif
#ifdef HAVE_SCALAPACK
    integer              :: ilbasis, jlbasis, proc(1:2), dest(1:2)
    integer              :: nevec_found, liwork, ii, node
    FLOAT,   allocatable :: gap(:)
    integer, allocatable :: iclustr(:)
    integer, allocatable :: send_count(:), send_disp(:), recv_count(:), recv_disp(:), recv_pos(:, :, :)
    R_TYPE,  allocatable :: send_buffer(:, :), recv_buffer(:, :)
    FLOAT                :: orfac
#ifdef R_TCOMPLEX
    FLOAT                :: rtmp(3) !< size must be at least 3 according to ScaLAPACK
    integer              :: lrwork
#endif
#endif
    type(profile_t), save :: prof

    PUSH_SUB(X(lcao_alt_wf).diagonalization)

    SAFE_ALLOCATE(ifail(1:this%nbasis))

    call profiling_in(prof, "LCAO_DIAG")

    if(this%parallel) then
#ifdef HAVE_SCALAPACK            
      SAFE_ALLOCATE(iclustr(1:2*st%dom_st_proc_grid%nprocs))
      SAFE_ALLOCATE(gap(1:st%dom_st_proc_grid%nprocs))

      ! This means that we do not want reorthogonalization of close
      ! eigenvectors
      orfac = M_ZERO

#ifdef R_TREAL
      call scalapack_sygvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = overlap(1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = orfac, &
        z = levec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = tmp(1), lwork = -1, iwork = liwork, liwork = -1, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)
#else
      call scalapack_hegvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = overlap(1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = orfac, &
        z = levec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = tmp(1), lwork = -1, rwork = rtmp(1), lrwork = -1, iwork = liwork, liwork = -1, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)
#endif

      if(info /= 0) then
        write(message(1), '(a,i4,a)') &
          'Workspace query for LCAO parallel diagonalization failed. ScaLAPACK returned info code ', info, '.'
        call messages_warning(1)
      end if

      lwork = nint(R_REAL(tmp(1)))

      SAFE_ALLOCATE(work(1:lwork))
      SAFE_ALLOCATE(iwork(1:liwork))

      SAFE_DEALLOCATE_A(ifail)
      SAFE_ALLOCATE(ifail(1:this%nbasis))

#ifdef R_TREAL
      call scalapack_sygvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = overlap(1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = orfac, &
        z = levec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = work(1), lwork = lwork, iwork = iwork(1), liwork = liwork, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)
#else
      lrwork = nint(rtmp(1))
      SAFE_ALLOCATE(rwork(1:lrwork))

      call scalapack_hegvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = overlap(1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = orfac, &
        z = levec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = work(1), lwork = lwork, rwork = rwork(1), lrwork = lrwork, iwork = iwork(1), liwork = liwork, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)
      SAFE_DEALLOCATE_A(rwork)
#endif

      if(info /= 0) then
        write(message(1), '(a,i4,a)') 'LCAO parallel diagonalization failed. ScaLAPACK returned info code ', info, '.'
        call messages_warning(1)
      end if

      ! Now we have to rearrange the data between processors. We have
      ! the data in levec, distributed according to ScaLAPACK and we
      ! need to copy it to evec, distributed by columns according to
      ! state parallelization.

      SAFE_ALLOCATE(send_count(1:st%dom_st_mpi_grp%size))
      SAFE_ALLOCATE(recv_count(1:st%dom_st_mpi_grp%size))

      ! First we count the number of points (to allocate the buffers)
      send_count = 0
      recv_count = 0
      do ibasis = 1, this%nbasis
        do jbasis = 1, st%nst

          call lcao_local_index(this, ibasis, jbasis, ilbasis, jlbasis, proc(1), proc(2))

          if(all(proc == this%myroc)) then
            ! we have to send the point

            ! we only count points per column (we have to send the
            ! same number to all points in a row)
            INCR(send_count(st%node(jbasis) + 1), 1)
          end if

          if(st%node(jbasis) == this%myroc(2)) then
            ! we have to receive
            call MPI_Cart_rank(st%dom_st_mpi_grp%comm, proc(1), node, mpi_err)
            INCR(recv_count(node + 1), 1)
          end if

        end do
      end do

      SAFE_ALLOCATE(send_buffer(1:max(1, maxval(send_count)), 1:st%dom_st_mpi_grp%size))
      SAFE_ALLOCATE(recv_pos(1:2, max(1, maxval(recv_count)), 1:st%dom_st_mpi_grp%size))

      send_count = 0
      recv_count = 0
      do ibasis = 1, this%nbasis
        do jbasis = 1, st%nst

          call lcao_local_index(this, ibasis, jbasis, ilbasis, jlbasis, proc(1), proc(2))

          if(all(proc == this%myroc)) then
            ! we have to send the point
            dest(2) = st%node(jbasis)

            do ii = 1, this%nproc(1)
              dest(1) = ii - 1

              ! get the node id from coordinates
              call MPI_Cart_rank(st%dom_st_mpi_grp%comm, dest(1), node, mpi_err)
              INCR(node, 1)
              INCR(send_count(node), 1)
              send_buffer(send_count(node), node) = levec(ilbasis, jlbasis)
            end do

          end if

          if(st%node(jbasis) == this%myroc(2)) then
            ! we have to receive
            call MPI_Cart_rank(st%dom_st_mpi_grp%comm, proc(1), node, mpi_err)
            INCR(node, 1)

            INCR(recv_count(node), 1)
            ! where do we put it once received?
            recv_pos(1, recv_count(node), node) = ibasis
            recv_pos(2, recv_count(node), node) = jbasis
          end if

        end do
      end do

      SAFE_ALLOCATE(recv_buffer(1:max(1, maxval(recv_count)), 1:st%dom_st_mpi_grp%size))

      SAFE_ALLOCATE(send_disp(1:st%dom_st_mpi_grp%size))
      SAFE_ALLOCATE(recv_disp(1:st%dom_st_mpi_grp%size))

      do node = 1, st%dom_st_mpi_grp%size
        send_disp(node) = ubound(send_buffer, dim = 1)*(node - 1)
        recv_disp(node) = ubound(recv_buffer, dim = 1)*(node - 1)
      end do

      call MPI_Alltoallv(send_buffer(1, 1), send_count(1), send_disp(1), R_MPITYPE, &
        recv_buffer(1, 1), recv_count(1), recv_disp(1), R_MPITYPE, &
        st%dom_st_mpi_grp%comm, mpi_err)

      do node = 1, st%dom_st_mpi_grp%size
        do ii = 1, recv_count(node)
          evec(recv_pos(1, ii, node), recv_pos(2, ii, node)) = recv_buffer(ii, node)
        end do
      end do

      SAFE_DEALLOCATE_A(iclustr)
      SAFE_DEALLOCATE_A(gap)

      SAFE_DEALLOCATE_A(send_disp)
      SAFE_DEALLOCATE_A(send_count)      
      SAFE_DEALLOCATE_A(send_buffer)
      SAFE_DEALLOCATE_A(recv_pos)
      SAFE_DEALLOCATE_A(recv_disp)
      SAFE_DEALLOCATE_A(recv_count)      
      SAFE_DEALLOCATE_A(recv_buffer)

#endif /* HAVE_SCALAPACK */
    else

      SAFE_ALLOCATE(iwork(1:5*this%nbasis))

#ifdef R_TREAL    
      call lapack_sygvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = tmp(1), lwork = -1, iwork = iwork(1), ifail = ifail(1), info = info)
#else

      SAFE_ALLOCATE(rwork(1:7*this%nbasis))

      call lapack_hegvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = tmp(1), lwork = -1, rwork = rwork(1), iwork = iwork(1), ifail = ifail(1), info = info)

#endif
      if(info /= 0) then
        write(message(1), '(a,i4,a)') 'Workspace query for LCAO diagonalization failed. LAPACK returned info code ', info, '.'
        call messages_warning(1)
      end if

      lwork = nint(R_REAL(tmp(1)))

      SAFE_ALLOCATE(work(1:lwork))

#ifdef R_TREAL
      call lapack_sygvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = work(1), lwork = lwork, iwork = iwork(1), ifail = ifail(1), info = info)
#else
      call lapack_hegvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = this%diag_tol, &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = work(1), lwork = lwork,  rwork = rwork(1), iwork = iwork(1), ifail = ifail(1), info = info)

      SAFE_DEALLOCATE_A(rwork)
#endif

      if(info /= 0) then
        write(message(1), '(a,i4,a)') 'LCAO diagonalization failed. LAPACK returned info code ', info, '.'
        call messages_warning(1)
      end if

#ifdef HAVE_MPI
    ! the eigenvectors are not unique due to phases and degenerate subspaces, but 
    ! they must be consistent among processors in domain parallelization
    if(gr%mesh%parallel_in_domains) &
      call MPI_Bcast(hamiltonian(1, 1), size(hamiltonian), R_MPITYPE, 0, gr%mesh%mpi_grp%comm, mpi_err)
#endif

    end if

    SAFE_DEALLOCATE_A(iwork)
    SAFE_DEALLOCATE_A(work)
    SAFE_DEALLOCATE_A(ifail)

    call profiling_out(prof)

    POP_SUB(X(lcao_alt_wf).diagonalization)
  end subroutine diagonalization

end subroutine X(lcao_alt_wf)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
