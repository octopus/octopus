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
! This routine fills state psi with an atomic orbital -- provided
! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(lcao_atomic_orbital) (this, iorb, m, hm, geo, sb, psi, spin_channel)
  type(lcao_t),             intent(in)    :: this
  integer,                  intent(in)    :: iorb
  type(mesh_t),             intent(in)    :: m
  type(simul_box_t),        intent(in)    :: sb
  type(hamiltonian_t),      intent(in)    :: hm
  type(geometry_t), target, intent(in)    :: geo
  R_TYPE,                   intent(inout) :: psi(:, :)
  integer,                  intent(in)    :: spin_channel

  type(species_t), pointer :: s
  type(periodic_copy_t)   :: pc
  integer :: icell, idim, iatom, jj, ip
  FLOAT :: pos(MAX_DIM)
  FLOAT, allocatable :: ao(:)
  type(profile_t), save :: prof

  call profiling_in(prof, "ATOMIC_ORBITAL")
  call push_sub('lcao_inc.Xlcao_atomic_orbital')

  ASSERT(iorb >= 1)
  ASSERT(iorb <= this%maxorbs)

  psi(1:m%np, 1:hm%d%dim) = R_TOTYPE(M_ZERO)

  iatom = this%atom(iorb)
  jj = this%level(iorb)
  idim = this%ddim(iorb)
  s => geo%atom(iatom)%spec
  ASSERT(jj <= species_niwfs(s))

  SAFE_ALLOCATE(ao(1:m%np))

  if (.not. simul_box_is_periodic(sb)) then

    call species_get_orbital(s, m, jj, calc_dim, max(spin_channel, idim), geo%atom(iatom)%x, ao)

    do ip = 1, m%np
      psi(ip, idim) = ao(ip)
    end do

  else

    call periodic_copy_init(pc, sb, geo%atom(iatom)%x, range = species_get_iwf_radius(s, jj, spin_channel))
    do icell = 1, periodic_copy_num(pc)
      pos = periodic_copy_position(pc, sb, icell)

      call species_get_orbital(s, m, jj, calc_dim, max(spin_channel, idim), pos, ao)
      
      do ip = 1, m%np
        psi(ip, idim) = psi(ip, idim) + ao(ip)
      end do

    end do
    call periodic_copy_end(pc)

  end if

  SAFE_DEALLOCATE_A(ao)

  call pop_sub()
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
  call push_sub('lcao_inc.Xlcao_wf')
  
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
  max = this%norbs*(this%norbs + 1)/2

  message(1) = "Info: Getting hamiltonian matrix elements"
  call write_info(1)

  if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max)

  do n1 = 1, this%norbs
    
    do ispin = 1, st%d%spin_channels
      call get_ao(n1, ispin, lcaopsi(:, :, ispin), use_psi = .true.)
    end do

    do ik = kstart, kend
      ispin = states_dim_get_spin_index(st%d, ik)
      call X(hamiltonian_apply)(hm, gr, lcaopsi(:, :, ispin), hpsi(:, :, ik), n1, ik)
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
    call MPI_Allgatherv(tmp(:, kstart:), st%nst*(kend - kstart + 1), MPI_FLOAT, &
         st%eigenval, st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
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

  call pop_sub()

contains 

  subroutine init_orbitals()
    integer :: iorb, ispin, ist, ik, size
    R_TYPE, allocatable :: ao(:, :)

    call push_sub('lcao_inc.Xlcao_wf.init_orbitals')

    ! We calculate the atomic orbitals first. To save memory we put
    ! all the orbitals we can in the part of st%Xpsi that we are going
    ! to overwrite and then the rest is stored in a single precision
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

      size = (this%norbs - iorb + 1)*st%d%spin_channels
      write(message(1), '(a, i5, a)') "Info: Extra storage for ", size, " orbitals will be allocated."
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

    call pop_sub()

  end subroutine init_orbitals

  subroutine get_ao(iorb, ispin, ao, use_psi)
    integer, intent(in)    :: iorb
    integer, intent(in)    :: ispin
    R_TYPE,  intent(out)   :: ao(:, :)
    logical, intent(in)    :: use_psi

    integer :: idim

    call push_sub('lcao_inc.Xlcao_wf.get_ao')

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

    call pop_sub()

  end subroutine get_ao

end subroutine X(lcao_wf)

! ---------------------------------------------------------
subroutine X(lcao_wf2) (this, st, gr, geo, hm, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: start

  integer :: iatom, jatom, ik, ist, idim, ip, ispin
  integer :: nbasis, ibasis, jbasis, iorb, jorb, maxorb, norbs
  R_TYPE, allocatable :: hamiltonian(:, :), overlap(:, :)
  R_TYPE, allocatable :: psii(:, :, :), hpsi(:, :, :)
  FLOAT, allocatable :: ev(:)
  FLOAT, allocatable :: radius(:)
  type(submesh_t), allocatable :: sphere(:)
  integer, allocatable :: basis_atom(:), basis_orb(:), atom_orb_basis(:, :)
  type(batch_t) :: hpsib, psib
  type(batch_t), allocatable :: orbitals(:)
  FLOAT :: dist2, lapdist, maxradius
  type(profile_t), save :: prof_orbitals, prof_matrix, prof_wavefunction

#ifdef HAVE_MPI
  R_TYPE, allocatable :: tmp(:, :)
  type(profile_t), save :: commprof, comm2prof
#endif

  maxorb = 0
  nbasis = 0
  do iatom = 1, geo%natoms
    maxorb = max(maxorb, species_niwfs(geo%atom(iatom)%spec) )
    nbasis = nbasis + species_niwfs(geo%atom(iatom)%spec)
  end do

  write(message(1),'(a,i6,a)') 'Info: Performing LCAO calculation with ', nbasis, ' orbitals.'
  call write_info(1)

  SAFE_ALLOCATE(hamiltonian(1:nbasis, 1:nbasis))
  SAFE_ALLOCATE(overlap(1:nbasis, 1:nbasis))
  SAFE_ALLOCATE(ev(1:nbasis))

  SAFE_ALLOCATE(psii(1:gr%mesh%np_part, 1:st%d%dim, maxorb))
  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim, maxorb))
  SAFE_ALLOCATE(radius(1:geo%natoms))
  SAFE_ALLOCATE(sphere(1:geo%natoms))
  SAFE_ALLOCATE(orbitals(1:geo%natoms))
  SAFE_ALLOCATE(basis_atom(1:nbasis))
  SAFE_ALLOCATE(basis_orb(1:nbasis))
  SAFE_ALLOCATE(atom_orb_basis(1:geo%natoms, 1:maxorb))

  ! this is the extra distance that the laplacian adds to the localization radius
  lapdist = maxval(abs(gr%mesh%idx%enlarge)*gr%mesh%h)

  call profiling_in(prof_orbitals, "LCAO_ORBITALS")

  do iatom = 1, geo%natoms
    norbs = species_niwfs(geo%atom(iatom)%spec)
    maxradius = M_ZERO
    do iorb = 1, norbs
      maxradius = max(maxradius, species_get_iwf_radius(geo%atom(iatom)%spec, iorb, is = 1))
    end do

    radius(iatom) = maxradius

    ! initialize the radial grid
    call submesh_init_sphere(sphere(iatom), gr%mesh%sb, gr%mesh, geo%atom(iatom)%x, maxradius)
    call batch_init(orbitals(iatom), 1, norbs)
    call X(batch_new)(orbitals(iatom), 1, norbs, sphere(iatom)%ns)
  end do

  do ispin = 1, st%d%spin_channels

    if(st%d%spin_channels > 1) then
      write(message(1), '(a,i1)') 'Info: LCAO for spin channel ', ispin
      call write_info(1)
    end if

    message(1) = "Info: Calculating atomic orbitals"
    call write_info(1)

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, nbasis)

    ibasis = 0
    do iatom = 1, geo%natoms
      norbs = species_niwfs(geo%atom(iatom)%spec)

      do iorb = 1, norbs
        ibasis = ibasis + 1
        atom_orb_basis(iatom, iorb) = ibasis
        basis_atom(ibasis) = iatom
        basis_orb(ibasis) = iorb
      end do

      do iorb = 1, norbs
        ! allocate and calculate the orbitals
        call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere(iatom), iorb, st%d%dim, ispin, &
          geo%atom(iatom)%x, orbitals(iatom)%states(iorb)%dpsi(:, 1))
      end do

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ibasis, nbasis)
    end do

    if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

    call profiling_out(prof_orbitals)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle

      message(1) = 'Info: Calculating matrix elements'
      if(st%d%nik > st%d%spin_channels) write(message(1), '(a,a,i5)') message(1), ' for k-point ', ik
      call write_info(1)

      call profiling_in(prof_matrix, "LCAO_MATRIX")

      hamiltonian = R_TOTYPE(M_ZERO)
      overlap = R_TOTYPE(M_ZERO)

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, geo%natoms)

      do iatom = 1, geo%natoms
        norbs = species_niwfs(geo%atom(iatom)%spec)

        psii = M_ZERO

        call batch_init( psib, st%d%dim, atom_orb_basis(iatom, 1), atom_orb_basis(iatom, norbs), psii)
        call batch_init(hpsib, st%d%dim, atom_orb_basis(iatom, 1), atom_orb_basis(iatom, norbs), hpsi)

        call X(submesh_batch_add)(sphere(iatom), orbitals(iatom), psib)
        call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)


        do jatom = 1, geo%natoms
          if(jatom < iatom) cycle


          dist2 = sum((geo%atom(iatom)%x(1:MAX_DIM) - geo%atom(jatom)%x(1:MAX_DIM))**2)

          if(dist2 > (radius(iatom) + radius(jatom) + lapdist)**2) cycle

          ibasis = atom_orb_basis(iatom, 1)
          jbasis = atom_orb_basis(jatom, 1)

          call X(submesh_batch_dotp_matrix)(sphere(jatom),  hpsib, orbitals(jatom), hamiltonian(ibasis:, jbasis:), reduce= .false.)

          if(dist2 > (radius(iatom) + radius(jatom))**2) cycle

          call X(submesh_batch_dotp_matrix)(sphere(jatom), psib, orbitals(jatom), overlap(ibasis:, jbasis:), reduce= .false.)

          ! the other half of the matrix
          do iorb = 1, norbs
            ibasis = atom_orb_basis(iatom, iorb)
            do jorb = 1, species_niwfs(geo%atom(jatom)%spec)
              jbasis = atom_orb_basis(jatom, jorb)
              hamiltonian(jbasis, ibasis) = R_CONJ(hamiltonian(ibasis, jbasis))
              overlap(jbasis, ibasis) = R_CONJ(overlap(ibasis, jbasis))
            end do
          end do

        end do

        call batch_end(psib)
        call batch_end(hpsib)

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(iatom, geo%natoms)
      end do

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

#ifdef HAVE_MPI

      if(gr%mesh%parallel_in_domains) then
        call profiling_in(comm2prof, "LCAO_REDUCE")
#ifndef HAVE_MPI2
        SAFE_ALLOCATE(tmp(1:nbasis, 1:nbasis))
        tmp = hamiltonian
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(tmp), hamiltonian, nbasis**2, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
#ifndef HAVE_MPI2
        tmp = overlap
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(tmp), overlap, nbasis**2, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(tmp)
        call profiling_out(comm2prof)
      end if
#endif

      call lalg_geneigensolve(nbasis, hamiltonian, overlap, ev)

      call profiling_out(prof_matrix)

      call profiling_in(prof_wavefunction, "LCAO_WAVEFUNCTIONS")

      message(1) = "Info: Generating wave-functions"
      if(st%d%nik > st%d%spin_channels) write(message(1), '(a,a,i5)') message(1), ' for k-point ', ik
      call write_info(1)

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, nbasis)


      do ist = st%st_start, st%st_end
        st%eigenval(ist, ik) = ev(ist)
        do idim = 1, st%d%dim
          forall(ip = 1:gr%mesh%np) st%X(psi)(ip, idim, ist, ik) = R_TOTYPE(M_ZERO)
        end do
      end do

      call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%X(psi)(:, :, :, ik))

      ibasis = 0
      do iatom = 1, geo%natoms
        norbs = species_niwfs(geo%atom(iatom)%spec)
        call X(submesh_batch_add_matrix)(sphere(iatom), hamiltonian(ibasis + 1:, st%st_start:), orbitals(iatom), psib)
        ibasis = ibasis + norbs
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ibasis, nbasis)
      end do

      call batch_end(psib)

      if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')
      call profiling_out(prof_wavefunction)
    end do
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
  SAFE_DEALLOCATE_A(ev)
  SAFE_DEALLOCATE_A(basis_atom)
  SAFE_DEALLOCATE_A(sphere)
  SAFE_DEALLOCATE_A(radius)

end subroutine X(lcao_wf2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
