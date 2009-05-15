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
  ASSERT(jj <= s%niwfs)

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

  integer :: iatom, jatom, ik, ist, idim, ip
  integer :: nbasis, ibasis, jbasis, iorb, jorb, maxorb
  R_TYPE, allocatable :: hamiltonian(:, :), overlap(:, :)
  R_TYPE, allocatable :: psii(:, :), hpsi(:, :)
  FLOAT, allocatable :: ev(:)
  FLOAT, pointer :: orb(:, :)
  FLOAT, allocatable :: radius(:)
  type(submesh_t), allocatable :: sphere(:)
  integer, allocatable :: basis_index(:, :)
  type(batch_t) :: orbitals
  FLOAT :: dist2, lapdist

#ifdef HAVE_MPI
  R_TYPE, allocatable :: tmp(:, :)
  type(profile_t), save :: commprof, comm2prof
#endif

  maxorb = 0
  nbasis = 0
  do iatom = 1, geo%natoms
    maxorb = max(maxorb, geo%atom(iatom)%spec%niwfs)
    nbasis = nbasis + geo%atom(iatom)%spec%niwfs
  end do
  
  SAFE_ALLOCATE(hamiltonian(1:nbasis, 1:nbasis))
  SAFE_ALLOCATE(overlap(1:nbasis, 1:nbasis))
  SAFE_ALLOCATE(ev(1:nbasis))
  
  SAFE_ALLOCATE(psii(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(radius(1:nbasis))
  SAFE_ALLOCATE(sphere(1:nbasis))
  SAFE_ALLOCATE(basis_index(1:geo%natoms, 1:maxorb))

  ! this is the extra distance that the laplacian adds to the localization radius
  lapdist = maxval(abs(gr%mesh%idx%enlarge)*gr%mesh%h)

  call batch_init(orbitals, 1, nbasis)

  ibasis = 0
  do iatom = 1, geo%natoms
    do iorb = 1, geo%atom(iatom)%spec%niwfs
      ibasis = ibasis + 1
      basis_index(iatom, iorb) = ibasis

      ! initialize the radial grid
      radius(ibasis) = species_get_iwf_radius(geo%atom(iatom)%spec, iorb, is = 1)
      call submesh_init_sphere(sphere(ibasis), gr%mesh%sb, gr%mesh, geo%atom(iatom)%x, radius(ibasis))

      ! allocate and calculate the orbitals
      call dbatch_new_state(orbitals, ibasis, sphere(ibasis)%ns)
      call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere(ibasis), iorb, st%d%dim, 1, &
           geo%atom(iatom)%x, orbitals%states(ibasis)%dpsi(:, 1))

    end do
  end do


  do ik = 1, st%d%nik
  
    hamiltonian = R_TOTYPE(M_ZERO)
    overlap = R_TOTYPE(M_ZERO)

    do iatom = geo%atoms%start, geo%atoms%end
      do iorb = 1, geo%atom(iatom)%spec%niwfs
        ibasis = basis_index(iatom, iorb)

        do idim = 1,st%d%dim 
          forall(ip = 1:gr%mesh%np) psii(ip, idim) = M_ZERO
          call submesh_add_to_mesh(sphere(ibasis), orbitals%states(ibasis)%dpsi(:, 1), psii(:, idim))
        end do

        call X(hamiltonian_apply)(hm, gr, psii, hpsi, ibasis, ik)

        do jatom = iatom, geo%natoms
          do jorb = 1, geo%atom(jatom)%spec%niwfs
            jbasis = basis_index(jatom, jorb)
            
            dist2 = sum((geo%atom(iatom)%x(1:MAX_DIM) - geo%atom(jatom)%x(1:MAX_DIM))**2)

            if(dist2 > (radius(ibasis) + radius(jbasis) + lapdist)**2) cycle

            hamiltonian(jbasis, ibasis) = &
                 submesh_to_mesh_dotp(sphere(jbasis), st%d%dim, orbitals%states(jbasis)%dpsi(:, 1), hpsi, reduce = .false.)
            hamiltonian(ibasis, jbasis) = R_CONJ(hamiltonian(jbasis, ibasis))

            if(dist2 > (radius(ibasis) + radius(jbasis))**2) cycle
            
            overlap(jbasis, ibasis) = &
                 submesh_to_mesh_dotp(sphere(jbasis), st%d%dim, orbitals%states(jbasis)%dpsi(:, 1), psii, reduce = .false.)
            overlap(ibasis, jbasis) = R_CONJ(overlap(jbasis, ibasis))

          end do
        end do
        
      end do
    end do

#ifdef HAVE_MPI
    if(geo%atoms%parallel) then
      call profiling_in(commprof, "LCAO_BCAST")

      !clearly this can be improved by using MPI_Allgatherv
      do iatom = 1, geo%natoms
        do iorb = 1, geo%atom(iatom)%spec%niwfs
          ibasis = basis_index(iatom, iorb)

          do jatom = iatom, geo%natoms
            do jorb = 1, geo%atom(jatom)%spec%niwfs
              jbasis = basis_index(jatom, jorb)

              call MPI_Bcast(hamiltonian(jbasis, ibasis), 1, R_MPITYPE, geo%atoms%node(iatom), geo%atoms%mpi_grp%comm, mpi_err)
              hamiltonian(ibasis, jbasis) = R_CONJ(hamiltonian(jbasis, ibasis))

              call MPI_Bcast(overlap(jbasis, ibasis), 1, R_MPITYPE, geo%atoms%node(iatom), geo%atoms%mpi_grp%comm, mpi_err)
              overlap(ibasis, jbasis) = R_CONJ(overlap(jbasis, ibasis))

            end do
          end do
        end do
      end do

      call profiling_out(commprof)
    end if

    if(gr%mesh%parallel_in_domains) then
      call profiling_in(comm2prof, "LCAO_REDUCE")
      SAFE_ALLOCATE(tmp(1:nbasis, 1:nbasis))
      tmp = hamiltonian
      call MPI_Allreduce(tmp, hamiltonian, nbasis**2, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
      tmp = overlap
      call MPI_Allreduce(tmp, overlap, nbasis**2, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
      SAFE_DEALLOCATE_A(tmp)
      call profiling_out(comm2prof)
    end if
#endif

    call lalg_geneigensolve(nbasis, hamiltonian, overlap, ev)

    do ist = st%st_start, st%st_end
      forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) st%X(psi)(ip, idim, ist, ik) = R_TOTYPE(M_ZERO)
    end do
    
    ibasis = 0
    do iatom = 1, geo%natoms
      do iorb = 1, geo%atom(iatom)%spec%niwfs
        ibasis = ibasis + 1

        if(ibasis <= st%nst) st%eigenval(ibasis, ik) = ev(ibasis)

        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            call submesh_add_to_mesh(sphere(ibasis), orbitals%states(ibasis)%dpsi(:, 1), &
                 st%X(psi)(:, idim, ist, ik), factor = hamiltonian(ibasis, ist))
          end do
        end do
        
      end do
    end do

  end do

  do ibasis = 1, nbasis
    call dbatch_delete_state(orbitals, ibasis)
    call submesh_end(sphere(ibasis))
  end do
  call batch_end(orbitals)

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(hpsi)  
  SAFE_DEALLOCATE_A(hamiltonian)
  SAFE_DEALLOCATE_A(overlap)
  SAFE_DEALLOCATE_A(ev)
  SAFE_DEALLOCATE_A(basis_index)
  SAFE_DEALLOCATE_A(sphere)
  SAFE_DEALLOCATE_A(radius)

end subroutine X(lcao_wf2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
