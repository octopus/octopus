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
subroutine X(lcao_atomic_orbital) (this, iorb, mesh, hm, geo, sb, psi, spin_channel)
  type(lcao_t),             intent(in)    :: this
  integer,                  intent(in)    :: iorb
  type(mesh_t),             intent(in)    :: mesh
  type(simul_box_t),        intent(in)    :: sb
  type(hamiltonian_t),      intent(in)    :: hm
  type(geometry_t), target, intent(in)    :: geo
  R_TYPE,                   intent(inout) :: psi(:, :)
  integer,                  intent(in)    :: spin_channel

  type(species_t), pointer :: spec
  type(periodic_copy_t) :: pc
  integer :: icell, idim, iatom, jj, ip
  FLOAT :: pos(MAX_DIM), radius
  FLOAT, allocatable :: ao(:)
  type(profile_t), save :: prof

  call profiling_in(prof, "ATOMIC_ORBITAL")
  PUSH_SUB(X(lcao_atomic_orbital))

  ASSERT(iorb >= 1)
  ASSERT(iorb <= this%maxorbs)

  psi(1:mesh%np, 1:hm%d%dim) = R_TOTYPE(M_ZERO)

  iatom = this%atom(iorb)
  jj = this%level(iorb)
  idim = this%ddim(iorb)
  spec => geo%atom(iatom)%spec
  ASSERT(jj <= species_niwfs(spec))

  SAFE_ALLOCATE(ao(1:mesh%np))

  if (.not. simul_box_is_periodic(sb)) then

    call species_get_orbital(spec, mesh, jj, mesh%sb%dim, max(spin_channel, idim), geo%atom(iatom)%x, ao)

    do ip = 1, mesh%np
      psi(ip, idim) = ao(ip)
    end do

  else

    radius = min(species_get_iwf_radius(spec, jj, spin_channel), maxval(sb%lsize))

    call periodic_copy_init(pc, sb, geo%atom(iatom)%x, range = radius)
    do icell = 1, periodic_copy_num(pc)
      pos = periodic_copy_position(pc, sb, icell)

      call species_get_orbital(spec, mesh, jj, mesh%sb%dim, max(spin_channel, idim), pos, ao)
      
      do ip = 1, mesh%np
        psi(ip, idim) = psi(ip, idim) + ao(ip)
      end do

    end do
    call periodic_copy_end(pc)

  end if

  SAFE_DEALLOCATE_A(ao)

  POP_SUB(X(lcao_atomic_orbital))
  call profiling_out(prof)

end subroutine X(lcao_atomic_orbital)


! ---------------------------------------------------------
subroutine X(lcao_wf) (this, st, gr, geo, hm, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: start

  integer :: nst, ik, n1, n2, idim, lcao_start, ie, max
  R_TYPE, allocatable :: hpsi(:, :, :), overlap(:, :, :)
  FLOAT, allocatable :: ev(:)
  R_TYPE, allocatable :: hamilt(:, :, :), lcaopsi(:, :, :), lcaopsi2(:, :)
  integer :: kstart, kend, ispin

  integer, allocatable :: cst(:, :), ck(:, :)
  R_SINGLE, allocatable :: buff(:, :, :, :)
#ifdef HAVE_MPI
  FLOAT, allocatable :: tmp(:, :)
#endif

  PUSH_SUB(X(lcao_wf))
  
  write(message(1),'(a,i6,a)') 'Info: Performing initial LCAO calculation with ', &
       this%norbs,' orbitals.'
  call write_info(1)
          

  nst = st%nst
  kstart = st%d%kpt%start
  kend = st%d%kpt%end

  lcao_start = start
  if(st%parallel_in_states .and. st%st_start > start) lcao_start = st%st_start

  ! Allocation of variables

  SAFE_ALLOCATE(lcaopsi(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%spin_channels))
  SAFE_ALLOCATE(lcaopsi2(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(hpsi(gr%mesh%np, st%d%dim, kstart:kend))
  SAFE_ALLOCATE(hamilt(this%norbs, this%norbs, kstart:kend))
  SAFE_ALLOCATE(overlap(this%norbs, this%norbs, st%d%spin_channels))

  call init_orbitals()

  ie = 0
  max = this%norbs * (this%norbs + 1)/2

  message(1) = "Info: Getting Hamiltonian matrix elements."
  call write_info(1)

  if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max)

  do n1 = 1, this%norbs
    
    do ispin = 1, st%d%spin_channels
      call get_ao(n1, ispin, lcaopsi(:, :, ispin), use_psi = .true.)
    end do

    do ik = kstart, kend
      ispin = states_dim_get_spin_index(st%d, ik)
      call X(hamiltonian_apply)(hm, gr%der, lcaopsi(:, :, ispin), hpsi(:, :, ik), n1, ik)
    end do

    do n2 = n1, this%norbs
      do ispin = 1, st%d%spin_channels

        call get_ao(n2, ispin, lcaopsi2, use_psi = .true.)

        overlap(n1, n2, ispin) = X(mf_dotp)(gr%mesh, st%d%dim, lcaopsi(:, :, ispin), lcaopsi2)
        overlap(n2, n1, ispin) = R_CONJ(overlap(n1, n2, ispin))
        do ik = kstart, kend
          if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle
          hamilt(n1, n2, ik) = X(mf_dotp)(gr%mesh, st%d%dim, hpsi(:, :, ik), lcaopsi2)
          hamilt(n2, n1, ik) = R_CONJ(hamilt(n1, n2, ik))
        end do
      end do
      
      ie = ie + 1
    end do

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ie, max)
  end do

  if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

  SAFE_DEALLOCATE_A(hpsi)

  SAFE_ALLOCATE(ev(1:this%norbs))

  do ik =  kstart, kend
    ispin = states_dim_get_spin_index(st%d, ik)
    call lalg_geneigensolve(this%norbs, hamilt(1:this%norbs, 1:this%norbs, ik), overlap(:, :, ispin), ev)

    st%eigenval(lcao_start:nst, ik) = ev(lcao_start:nst)

    st%X(psi)(1:gr%mesh%np, 1:st%d%dim, lcao_start:st%st_end, ik) = R_TOTYPE(M_ZERO)
  end do
  
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

  ! Change of basis
  do n2 = 1, this%norbs
    do ispin = 1, st%d%spin_channels
      
      call get_ao(n2, ispin, lcaopsi2, use_psi = .false.)

      do ik =  kstart, kend
        if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle
        do idim = 1, st%d%dim
          do n1 = lcao_start, st%st_end
            call lalg_axpy(gr%mesh%np, hamilt(n2, n1, ik), lcaopsi2(:, idim), st%X(psi)(:, idim, n1, ik))
          end do
        end do
      end do
      
    end do
  end do

  SAFE_DEALLOCATE_A(ev)
  SAFE_DEALLOCATE_A(hamilt)
  SAFE_DEALLOCATE_A(overlap)
  
  SAFE_DEALLOCATE_A(buff)
  SAFE_DEALLOCATE_A(cst)
  SAFE_DEALLOCATE_A(ck)

  SAFE_DEALLOCATE_A(lcaopsi)
  SAFE_DEALLOCATE_A(lcaopsi2)

  POP_SUB(X(lcao_wf))

contains 


! ---------------------------------------------------------
  subroutine init_orbitals()
    integer :: iorb, ispin, ist, ik, size
    R_TYPE, allocatable :: ao(:, :)

    PUSH_SUB(X(lcao_wf).init_orbitals)

    ! We calculate the atomic orbitals first. To save memory we put
    ! all the orbitals we can in the part of st%Xpsi that we are going
    ! to overwrite and then the rest is stored in a single-precision
    ! buffer.

    SAFE_ALLOCATE(cst(1:this%norbs, 1:st%d%spin_channels))
    SAFE_ALLOCATE( ck(1:this%norbs, 1:st%d%spin_channels))

    ck = 0

    iorb = 1
    ispin = 1

    ! first store in st%Xpsi
    ist_loop: do ist = lcao_start, st%st_end
      do ik = kstart, kend

        cst(iorb, ispin) = ist
        ck(iorb, ispin) = ik

        call X(lcao_atomic_orbital)(this, iorb, gr%mesh, hm, geo, gr%sb, st%X(psi)(:, :, ist, ik), ispin)

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
      call write_info(1)

      SAFE_ALLOCATE(buff(1:gr%mesh%np, 1:st%d%dim, iorb:this%norbs, 1:st%d%spin_channels))
      SAFE_ALLOCATE(  ao(1:gr%mesh%np, 1:st%d%dim))
      
      do iorb = iorb, this%norbs
        do ispin = 1, st%d%spin_channels
          call X(lcao_atomic_orbital)(this, iorb, gr%mesh, hm, geo, gr%sb, ao, ispin)
          buff(1:gr%mesh%np, 1:st%d%dim, iorb, ispin) = ao(1:gr%mesh%np, 1:st%d%dim)
        end do
      end do
      
      SAFE_DEALLOCATE_A(ao)
    end if

    POP_SUB(X(lcao_wf).init_orbitals)

  end subroutine init_orbitals


! ---------------------------------------------------------
  subroutine get_ao(iorb, ispin, ao, use_psi)
    integer, intent(in)    :: iorb
    integer, intent(in)    :: ispin
    R_TYPE,  intent(out)   :: ao(:, :)
    logical, intent(in)    :: use_psi

    integer :: idim

    PUSH_SUB(X(lcao_wf).get_ao)

    if(ck(iorb, ispin) == 0) then
      ao(1:gr%mesh%np, 1:st%d%dim) = buff(1:gr%mesh%np, 1:st%d%dim, iorb, ispin)
    else
      if(use_psi) then
        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, cst(iorb, ispin), ck(iorb, ispin)), ao(:, idim))
        end do
      else
        call X(lcao_atomic_orbital)(this, iorb, gr%mesh, hm, geo, gr%sb, ao, ispin)
      end if
    end if

    POP_SUB(X(lcao_wf).get_ao)

  end subroutine get_ao

end subroutine X(lcao_wf)


! ---------------------------------------------------------
! The alternative implementation.
subroutine X(lcao_wf2) (this, st, gr, geo, hm, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: start

  integer :: iatom, jatom, ik, ist, idim, ip, ispin, nev
  integer :: ibasis, jbasis, iorb, jorb, norbs, dof
  R_TYPE, allocatable :: hamiltonian(:, :), overlap(:, :), aa(:, :), bb(:, :)
  R_TYPE, allocatable :: lhamiltonian(:, :, :), loverlap(:, :, :)
#ifdef HAVE_SCALAPACK
  integer :: prow, pcol, ilbasis, jlbasis
#endif
  R_TYPE, allocatable :: psii(:, :, :), hpsi(:, :, :)
  type(submesh_t), allocatable :: sphere(:)
  type(batch_t) :: hpsib, psib
  type(batch_t), allocatable :: orbitals(:)
  FLOAT, allocatable :: eval(:)
  R_TYPE, allocatable :: evec(:, :)
  FLOAT :: dist2
  type(profile_t), save :: prof_orbitals, prof_matrix, prof_wavefunction

  PUSH_SUB(X(lcao_wf2))

  write(message(1), '(a,i6,a)') 'Info: Performing LCAO calculation with ', this%nbasis, ' orbitals.'
  write(message(2), '(a)') ' '
  call write_info(2)

  if (this%nbasis < st%nst) then
    write(message(1), '(a)') 'Not enough atomic orbitals to initialize all states,'
    write(message(2), '(i6,a)') st%nst - this%nbasis, ' states will be randomized.'
    if(this%derivative) then
      call write_warning(2)
    else
      write(message(3), '(a)') 'You can duplicate the number of atomic orbitals by setting'
      write(message(4), '(a)') 'LCAOExtraOrbitals to yes.'
      call write_warning(4)
    end if
  end if

  if(.not. this%parallel) then
    SAFE_ALLOCATE(hamiltonian(1:this%nbasis, 1:this%nbasis))
    SAFE_ALLOCATE(overlap(1:this%nbasis, 1:this%nbasis))
  else
#ifdef HAVE_SCALAPACK
    SAFE_ALLOCATE(lhamiltonian(1:this%lsize(1), 1:this%lsize(2), 1:this%nproc(1)))
    SAFE_ALLOCATE(loverlap(1:this%lsize(1), 1:this%lsize(2), 1:this%nproc(1)))
#endif
  end if
  SAFE_ALLOCATE(aa(1:this%maxorb, 1:this%maxorb))
  SAFE_ALLOCATE(bb(1:this%maxorb, 1:this%maxorb))  
  SAFE_ALLOCATE(psii(1:gr%mesh%np_part, 1:st%d%dim, this%maxorb))
  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim, this%maxorb))
  SAFE_ALLOCATE(sphere(1:geo%natoms))
  SAFE_ALLOCATE(orbitals(1:geo%natoms))

  call profiling_in(prof_orbitals, "LCAO_ORBITALS")

  dof = 0
  do iatom = 1, geo%natoms
    norbs = species_niwfs(geo%atom(iatom)%spec)

    ! initialize the radial grid
    call submesh_init_sphere(sphere(iatom), gr%mesh%sb, gr%mesh, geo%atom(iatom)%x, this%radius(iatom))
    INCR(dof, sphere(iatom)%ns*this%mult*norbs)
    call batch_init(orbitals(iatom), 1, this%mult*norbs)
    call dbatch_new(orbitals(iatom), 1, this%mult*norbs, sphere(iatom)%ns)
  end do

  write(message(1), '(a,i8,a,i7,a)')'Info: LCAO requires', dof, ' values (', dof*8/1024**2, ' Mb) for atomic orbitals.'
  call write_info(1)

  do ispin = 1, st%d%spin_channels

    if(st%d%spin_channels > 1) then
      write(message(1), '(a)') ' '
      write(message(2), '(a,i1)') 'LCAO for spin channel ', ispin
      write(message(3), '(a)') ' '
      call write_info(3)
    end if

    message(1) = 'Generating atomic orbitals.'
    call write_info(1)

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, this%nbasis)

    ibasis = 0
    do iatom = 1, geo%natoms
      norbs = species_niwfs(geo%atom(iatom)%spec)

      do iorb = 1, norbs
        ! allocate and calculate the orbitals
        call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere(iatom), iorb, st%d%dim, ispin, &
          geo%atom(iatom)%x, orbitals(iatom)%states(iorb)%dpsi(:, 1))
        if(this%derivative) then
          call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere(iatom), iorb, st%d%dim, ispin, &
            geo%atom(iatom)%x, orbitals(iatom)%states(iorb + norbs)%dpsi(:, 1), derivative = .true.)
        end if
      end do

      INCR(ibasis, norbs*this%mult)
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ibasis, this%nbasis)
    end do

    if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

    call profiling_out(prof_orbitals)
    
    ! iterate over the kpoints for this spin
    do ik = st%d%kpt%start, st%d%kpt%end
      if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle

      if(st%d%nik > st%d%spin_channels) then
        write(message(1), '(a)') ' '
        write(message(2), '(a,i5)') 'LCAO for k-point ', states_dim_get_kpoint_index(st%d, ik)
        write(message(3), '(a)') ' '
        call write_info(3)
      end if

      message(1) = 'Calculating matrix elements.'
      call write_info(1)

      call profiling_in(prof_matrix, "LCAO_MATRIX")

      if(.not. this%parallel) then
        hamiltonian = R_TOTYPE(M_ZERO)
        overlap = R_TOTYPE(M_ZERO)
      else
        lhamiltonian = R_TOTYPE(M_ZERO)
        loverlap = R_TOTYPE(M_ZERO)
      end if

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, geo%natoms)

      do iatom = 1, geo%natoms
        norbs = this%mult*species_niwfs(geo%atom(iatom)%spec)

        psii = M_ZERO

        call batch_init( psib, st%d%dim, this%atom_orb_basis(iatom, 1), this%atom_orb_basis(iatom, norbs), psii)
        call batch_init(hpsib, st%d%dim, this%atom_orb_basis(iatom, 1), this%atom_orb_basis(iatom, norbs), hpsi)

        call X(submesh_batch_add)(sphere(iatom), orbitals(iatom), psib)
        call X(hamiltonian_apply_batch)(hm, gr%der, psib, hpsib, ik)

        do jatom = 1, geo%natoms
          ! we only calculate the upper triangle
          if(jatom < iatom) cycle

          dist2 = sum((geo%atom(iatom)%x(1:MAX_DIM) - geo%atom(jatom)%x(1:MAX_DIM))**2)

          if(dist2 > (this%radius(iatom) + this%radius(jatom) + this%lapdist)**2) cycle

          ibasis = this%atom_orb_basis(iatom, 1)
          jbasis = this%atom_orb_basis(jatom, 1)

          call X(submesh_batch_dotp_matrix)(sphere(jatom), hpsib, orbitals(jatom), aa, reduce = .false.)

          if(dist2 > (this%radius(iatom) + this%radius(jatom))**2) then 
            aa = M_ZERO
          else
            call X(submesh_batch_dotp_matrix)(sphere(jatom), psib, orbitals(jatom), bb, reduce = .false.)
          end if

          if(.not. this%parallel) then
            forall(iorb = 1:norbs, jorb = 1:norbs)
              hamiltonian(ibasis - 1 + iorb, jbasis - 1 + jorb) = aa(iorb, jorb)
              overlap(ibasis - 1 + iorb, jbasis - 1 + jorb) = bb(iorb, jorb)
            end forall
          else
#ifdef HAVE_SCALAPACK
            do iorb = 1, norbs
              do jorb = 1, norbs
                call lcao_local_index(this, ibasis - 1 + iorb,  jbasis - 1 + jorb, &
                  ilbasis, jlbasis, prow, pcol)
                if(pcol == this%myroc(2)) then
                  lhamiltonian(ilbasis, jlbasis, prow + 1) = aa(iorb, jorb)
                  loverlap(ilbasis, jlbasis, prow + 1) = bb(iorb, jorb)
                end if
              end do
            end do
#endif
          end if

        end do

        call batch_end(psib)
        call batch_end(hpsib)

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(iatom, geo%natoms)
      end do

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

      call reduce()

      message(1) = 'Diagonalizing Hamiltonian.'
      call write_info(1)

      call profiling_out(prof_matrix)

      ! the number of eigenvectors we need
      nev = min(this%nbasis, st%nst) 

      SAFE_ALLOCATE(eval(1:this%nbasis))
      SAFE_ALLOCATE(evec(1:this%nbasis, 1:this%nbasis))

      call diagonalization()

      call profiling_in(prof_wavefunction, "LCAO_WAVEFUNCTIONS")

      message(1) = 'Generating wavefunctions.'
      call write_info(1)

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, this%nbasis)

      ! set the eigenvalues
      st%eigenval(1:nev, ik) = eval(1:nev)
      st%eigenval(this%nbasis + 1:st%nst, ik) = HUGE(M_ONE)

      do ist = st%st_start, st%st_end
        if(ist <= this%nbasis) then
          do idim = 1, st%d%dim
            forall(ip = 1:gr%mesh%np) st%X(psi)(ip, idim, ist, ik) = R_TOTYPE(M_ZERO)
          end do
        else
          !we do not have enough basis functions, so randomize the missing orbitals
          call X(mf_random)(gr%mesh, st%X(psi)(:, 1, ist, ik))
          st%X(psi)(1:gr%mesh%np, 2:st%d%dim, ist, ik) = R_TOTYPE(M_ZERO)
        end if
      end do

      call batch_init(psib, st%d%dim, st%st_start,  min(st%st_end, this%nbasis), st%X(psi)(:, :, :, ik))

      ibasis = 1
      do iatom = 1, geo%natoms
        norbs = this%mult*species_niwfs(geo%atom(iatom)%spec)
        call X(submesh_batch_add_matrix)(sphere(iatom), evec(ibasis:, st%st_start:), orbitals(iatom), psib)
        ibasis = ibasis + norbs
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ibasis, this%nbasis)
      end do

      call batch_end(psib)

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')
      call profiling_out(prof_wavefunction)
    end do

    SAFE_DEALLOCATE_A(eval)
    SAFE_DEALLOCATE_A(evec)

  end do

  do iatom = 1, geo%natoms
    call submesh_end(sphere(iatom))
    call X(batch_delete)(orbitals(iatom))
    call batch_end(orbitals(iatom))
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(hpsi)  
  SAFE_DEALLOCATE_A(hamiltonian)
  SAFE_DEALLOCATE_A(overlap)
  SAFE_DEALLOCATE_A(sphere)

  POP_SUB(X(lcao_wf2))

contains

  subroutine diagonalization()
    
    integer :: neval_found
    integer :: info, lwork
    FLOAT   :: tmp
    R_TYPE, allocatable :: work(:)
    integer, allocatable :: iwork(:), ifail(:)
#ifdef HAVE_SCALAPACK
    integer :: nevec_found, liwork
    FLOAT, allocatable :: gap(:)
    integer, allocatable :: iclustr(:)
#endif

    SAFE_ALLOCATE(ifail(1:this%nbasis))
    
#ifdef R_TREAL

    if(this%parallel) then
#ifdef HAVE_SCALAPACK            
      SAFE_ALLOCATE(iclustr(1:2*st%dom_st_proc_grid%nprocs))
      SAFE_ALLOCATE(gap(1:st%dom_st_proc_grid%nprocs))

      call scalapack_sygvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = lhamiltonian(1, 1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = loverlap(1, 1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = CNST(1e-15), &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = -M_ONE, &
        z = evec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = tmp, lwork = -1, iwork = liwork, liwork = -1, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)

      ASSERT(info == 0)

      lwork = nint(tmp)

      SAFE_ALLOCATE(work(1:lwork))
      SAFE_ALLOCATE(iwork(1:liwork))

      call scalapack_sygvx(ibtype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = lhamiltonian(1, 1, 1), ia = 1, ja = 1, desca = this%desc(1), &
        b = loverlap(1, 1, 1), ib = 1, jb = 1, descb = this%desc(1), &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = CNST(1e-15), &
        m = neval_found, nz = nevec_found, w = eval(1), orfac = -M_ONE, &
        z = evec(1, 1), iz = 1, jz = 1, descz = this%desc(1), &
        work = work(1), lwork = lwork, iwork = iwork(1), liwork = liwork, &
        ifail = ifail(1), iclustr = iclustr(1), gap = gap(1), info = info)
      
      ASSERT(info == 0)

      call MPI_Bcast(evec(1, 1), this%nbasis**2, R_MPITYPE, 0, st%dom_st_mpi_grp%comm, mpi_err)

      SAFE_DEALLOCATE_A(iclustr)
      SAFE_DEALLOCATE_A(gap)
#endif
    else

      SAFE_ALLOCATE(iwork(1:5*this%nbasis))
    
      call lapack_sygvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = CNST(1e-15), &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = tmp, lwork = -1, iwork = iwork(1), ifail = ifail(1), info = info)

      ASSERT(info == 0)
      
      lwork = nint(tmp)
      
      SAFE_ALLOCATE(work(1:lwork))
      
      call lapack_sygvx(itype = 1, jobz = 'V', range = 'I', uplo = 'U', &
        n = this%nbasis, a = hamiltonian(1, 1), lda = this%nbasis, b = overlap(1, 1), ldb = this%nbasis, &
        vl = M_ZERO, vu = M_ONE, il = 1, iu = nev, abstol = CNST(1e-15), &
        m = neval_found, w = eval(1), z = evec(1, 1), ldz = this%nbasis, &
        work = work(1), lwork = lwork, iwork = iwork(1), ifail = ifail(1), info = info)
      
    end if
    
#else
    call lalg_geneigensolve(this%nbasis, hamiltonian, overlap, eval)
    evec(this%nbasis, 1:nev) = hamiltonian(this%nbasis, 1:nev)
#endif
    SAFE_DEALLOCATE_A(iwork)
    SAFE_DEALLOCATE_A(work)
    SAFE_DEALLOCATE_A(ifail)

  end subroutine diagonalization

  ! -----------------------------------------------

  subroutine reduce()
#ifdef HAVE_MPI
#ifdef HAVE_SCALAPACK
    R_TYPE, allocatable :: ltmp(:, :, :)
#endif
    R_TYPE, allocatable :: tmp(:, :)
    type(profile_t), save :: comm2prof
    integer :: size

    if(gr%mesh%parallel_in_domains) then
      call profiling_in(comm2prof, "LCAO_REDUCE")
      if(.not. this%parallel) then
        SAFE_ALLOCATE(tmp(1:this%nbasis, 1:this%nbasis))
        tmp = hamiltonian
        call MPI_Allreduce(tmp(1, 1), hamiltonian(1, 1), this%nbasis**2, R_MPITYPE, MPI_SUM, &
          gr%mesh%mpi_grp%comm, mpi_err)
        tmp = overlap
        call MPI_Allreduce(tmp(1, 1), overlap(1, 1), this%nbasis**2, R_MPITYPE, MPI_SUM, &
          gr%mesh%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(tmp)
      else

#ifdef HAVE_SCALAPACK
        SAFE_ALLOCATE(ltmp(1:this%lsize(1), 1:this%lsize(2), 1:this%nproc(1)))

        ! This does not work yet. Size is not the same in all nodes...
        size = this%lsize(1)*this%lsize(2)*this%nproc(1)

        call MPI_Barrier(st%dom_st_mpi_grp%comm, mpi_err)

        ltmp = lhamiltonian
        call MPI_Allreduce(ltmp(1, 1, 1), lhamiltonian(1, 1, 1), size, R_MPITYPE, MPI_SUM, &
          gr%mesh%mpi_grp%comm, mpi_err)
        ltmp = loverlap
        call MPI_Allreduce(ltmp(1, 1, 1), loverlap(1, 1, 1), size, R_MPITYPE, MPI_SUM, &
          gr%mesh%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(ltmp)
#endif
      end if
      call profiling_out(comm2prof)
    end if
#endif
  end subroutine reduce

end subroutine X(lcao_wf2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
