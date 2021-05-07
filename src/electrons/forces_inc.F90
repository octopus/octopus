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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

subroutine X(forces_gather)(ions, force)
  type(ions_t),     intent(in)    :: ions
  R_TYPE,           intent(inout) :: force(:, :)
  
  R_TYPE,  allocatable :: force_local(:, :)
  integer, allocatable :: recv_count(:), recv_displ(:)

  PUSH_SUB(X(forces_gather))

  call profiling_in(prof_comm, TOSTRING(X(FORCES_COMM)))
    
  ! each node has a piece of the force array, they have to be
  ! collected by all nodes
  
  ASSERT(ubound(force, dim = 1) == ions%space%dim)

  SAFE_ALLOCATE(recv_count(1:ions%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(recv_displ(1:ions%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(force_local(1:ions%space%dim, 1:max(1, ions%atoms_dist%nlocal)))
  
  recv_count(1:ions%atoms_dist%mpi_grp%size) = ions%space%dim*ions%atoms_dist%num(0:ions%atoms_dist%mpi_grp%size - 1)
  recv_displ(1:ions%atoms_dist%mpi_grp%size) = ions%space%dim*(ions%atoms_dist%range(1, 0:ions%atoms_dist%mpi_grp%size - 1) - 1)
  
  if(ions%atoms_dist%nlocal > 0) then
    force_local(1:ions%space%dim, 1:ions%atoms_dist%nlocal) = force(1:ions%space%dim, ions%atoms_dist%start:ions%atoms_dist%end)
  end if

#ifdef HAVE_MPI  
  call MPI_Allgatherv(&
    force_local(1, 1), ions%space%dim*ions%atoms_dist%nlocal, R_MPITYPE, &
    force(1, 1), recv_count(1), recv_displ(1), R_MPITYPE, &
    ions%atoms_dist%mpi_grp%comm, mpi_err)
#endif
  
  SAFE_DEALLOCATE_A(recv_count)
  SAFE_DEALLOCATE_A(recv_displ)
  SAFE_DEALLOCATE_A(force_local)

  call profiling_out(prof_comm)
  POP_SUB(X(forces_gather))
end subroutine X(forces_gather)

!---------------------------------------------------------------------------
subroutine X(forces_from_local_potential)(mesh, namespace, ions, ep, gdensity, force)
  type(mesh_t),                   intent(in)    :: mesh
  type(namespace_t),              intent(in)    :: namespace
  type(ions_t),                   intent(in)    :: ions
  type(epot_t),                   intent(in)    :: ep
  R_TYPE,                         intent(in)    :: gdensity(:, :)
  R_TYPE,                         intent(inout) :: force(:, :)

  FLOAT,  allocatable :: vloc(:)
  integer             :: ip, idir, iatom
  R_TYPE, allocatable  :: zvloc(:), force_tmp(:,:)
  type(profile_t), save :: prof
 
  PUSH_SUB(X(forces_from_local_potential))

  call profiling_in(prof, TOSTRING(X(FORCES_LOCAL_POT)))

  SAFE_ALLOCATE(vloc(1:mesh%np))
  SAFE_ALLOCATE(zvloc(1:mesh%np))

  SAFE_ALLOCATE(force_tmp(1:ions%space%dim, 1:ions%natoms))
  force_tmp = M_ZERO
  
  do iatom = ions%atoms_dist%start, ions%atoms_dist%end

    vloc(1:mesh%np) = M_ZERO
    call epot_local_potential(ep, namespace, ions%space, ions%latt, mesh, ions%atom(iatom)%species, &
      ions%pos(:, iatom), iatom, vloc)

    do ip = 1, mesh%np
      zvloc(ip) = vloc(ip)
    end do

    do idir = 1, ions%space%dim
      force_tmp(idir, iatom) = -X(mf_dotp)(mesh, zvloc, gdensity(:, idir), reduce = .false.)
    end do

  end do

  if(ions%atoms_dist%parallel) call X(forces_gather)(ions, force_tmp)
  !if(ions%atoms_dist%parallel .and. ions%atoms_dist%nlocal > 0) call X(forces_gather)(ions, force)

  if(mesh%parallel_in_domains) call mesh%allreduce(force_tmp)

  force(1:ions%space%dim, 1:ions%natoms) = force(1:ions%space%dim, 1:ions%natoms) + force_tmp(1:ions%space%dim, 1:ions%natoms)

  SAFE_DEALLOCATE_A(vloc)
  SAFE_DEALLOCATE_A(zvloc)
  SAFE_DEALLOCATE_A(force_tmp)

  call profiling_out(prof)

  POP_SUB(X(forces_from_local_potential))
end subroutine X(forces_from_local_potential)

!---------------------------------------------------------------------------
!> Ref: Kikuji Hirose, Tomoya Ono, Yoshitaka Fujimoto, and Shigeru Tsukamoto,
!! First-principles calculations in real-space formalism: Electronic configurations
!! and transport properties of nanostructures, Imperial College Press (2005)
!! Section 1.6, page 12
subroutine X(forces_from_potential)(gr, namespace, space, ions, hm, st, force, force_loc, force_nl, force_u)
  type(grid_t),                   intent(in)    :: gr
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(ions_t),                   intent(in)    :: ions
  type(hamiltonian_elec_t),       intent(in)    :: hm
  type(states_elec_t),            intent(in)    :: st
  FLOAT,                          intent(out)   :: force(:, :)
  FLOAT,                          intent(out)   :: force_loc(:, :)
  FLOAT,                          intent(out)   :: force_nl(:, :)
  FLOAT,                          intent(out)   :: force_u(:, :)


  type(symmetrizer_t) :: symmetrizer
  integer :: iatom, ist, iq, idim, idir, np, np_part, ikpoint, iop, ii, iatom_symm
  integer :: ib, maxst, minst, ip
  FLOAT :: ratom(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  FLOAT,  allocatable :: grad_rho(:, :)
  FLOAT,  allocatable :: force_psi(:)
  FLOAT, allocatable :: symmtmp(:, :)
  type(wfs_elec_t) :: psib, grad_psib(1:MAX_DIM)
  FLOAT :: kweight
  type(profile_t), save :: prof

  PUSH_SUB(X(forces_from_potential))

  call profiling_in(prof, TOSTRING(X(FORCES_FROM_POTENTIALS)))

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_psi(1:np, 1:space%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:space%dim))
  grad_rho = M_ZERO

  SAFE_ALLOCATE(force_psi(1:space%dim))

  force = M_ZERO
  force_loc = M_ZERO
  force_nl = M_ZERO
  force_u = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end

    ikpoint = st%d%get_kpoint_index(iq)
    if(st%d%kweights(iq) <= M_EPSILON) cycle

    do ib = st%group%block_start, st%group%block_end
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      call st%group%psib(ib, iq)%copy_to(psib, copy_data = .true.)
      if(hamiltonian_elec_apply_packed(hm)) call psib%do_pack()

      ! set the boundary conditions
      call boundaries_set(gr%der%boundaries, psib)

      ! set the phase for periodic systems
      if (allocated(hm%hm_base%phase)) then
        call hamiltonian_elec_base_phase(hm%hm_base, gr%mesh, gr%mesh%np_part, .false., psib)
      end if

      ! calculate the gradient
      do idir = 1, space%dim
        call psib%copy_to(grad_psib(idir))
        call X(derivatives_batch_perform)(gr%der%grad(idir), gr%der, psib, grad_psib(idir), set_bc = .false.)
      end do

      ! calculate the contribution to the density gradient
      call X(density_accumulate_grad)(gr, st, psib, grad_psib, grad_rho)

      ! the non-local potential contribution
      if(hm%hm_base%apply_projector_matrices .and. .not. accel_is_enabled() .and. &
        .not. (st%symmetrize_density .and. hm%kpoints%use_symmetries)) then

        call X(hamiltonian_elec_base_nlocal_force)(hm%hm_base, gr%mesh, st, gr%der%boundaries, iq, &
          space%dim, psib, grad_psib, force_nl)

      else 

        do ist = minst, maxst

          if(abs(st%occ(ist, iq)) <= M_EPSILON) cycle

          ! get the state and its gradient out of the batches (for the moment)
          do idim = 1, st%d%dim
            call batch_get_state(psib, (/ist, idim/), gr%mesh%np_part, psi(:, idim))
            do idir = 1, space%dim
              call batch_get_state(grad_psib(idir), (/ist, idim/), gr%mesh%np, grad_psi(:, idir, idim))
            end do
          end do

          call profiling_count_operations(np*st%d%dim*space%dim*(2 + R_MUL))

          if(st%symmetrize_density .and. hm%kpoints%use_symmetries) then

            ! We use that
            !
            ! \int dr f(Rr) V_iatom(r) \nabla f(R(v)) = R\int dr f(r) V_iatom(R*r) f(r)
            !
            ! and that the operator R should map the position of atom
            ! iatom to the position of some other atom iatom_symm, so that
            !
            ! V_iatom(R*r) = V_iatom_symm(r)
            !

            !We apply N symmetries, so we have to use the proper weight
            kweight = st%d%kweights(iq) / kpoints_get_num_symmetry_ops(hm%kpoints, ikpoint)

            do ii = 1, kpoints_get_num_symmetry_ops(hm%kpoints, ikpoint)

              iop = abs(kpoints_get_symmetry_ops(hm%kpoints, ikpoint, ii))
              !if(iop < 0 ) cycle !Time reversal symmetry

              do iatom = 1, ions%natoms
                if(projector_is_null(hm%ep%proj(iatom))) cycle

                !We find the atom that correspond to this one, once symmetry is applied
                ratom(1:space%dim) = symm_op_apply_inv_cart(gr%symm%ops(iop), ions%pos(:, iatom))

                ratom(1:ions%space%dim) = ions%latt%fold_into_cell(ratom(1:ions%space%dim))

                ! find iatom_symm
                do iatom_symm = 1, ions%natoms
                  if(all(abs(ratom(1:space%dim) - ions%pos(:, iatom_symm)) < CNST(1.0e-5))) exit
                end do

                if(iatom_symm > ions%natoms) then
                  write(message(1),'(a,i6)') 'Internal error: could not find symmetric partner for atom number', iatom
                  write(message(2),'(a,i3,a)') 'with symmetry operation number ', iop, '.'
                  call messages_fatal(2, namespace=namespace)
                end if

                do idir = 1, space%dim
                  force_psi(idir) = - M_TWO * kweight * st%occ(ist, iq) * &
      R_REAL(X(projector_matrix_element)(hm%ep%proj(iatom_symm), gr%der%boundaries, st%d%dim, iq, psi, grad_psi(:, idir, :)))
                end do

                ! We convert the force to Cartesian coordinates before symmetrization
                ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
                if (space%is_periodic() .and. gr%sb%latt%nonorthogonal ) then 
                  force_psi(1:space%dim) = matmul(gr%sb%latt%klattice_primitive(1:space%dim, 1:space%dim), force_psi(1:space%dim))
                end if

                !Let us now apply the symmetry to the force
                !Note: here we are working with reduced quantities
                force_nl(1:space%dim, iatom) = force_nl(1:space%dim, iatom) + symm_op_apply_cart(gr%symm%ops(iop), force_psi)

              end do

            end do

          else

            ! iterate over the projectors
            do iatom = 1, ions%natoms
              if(projector_is_null(hm%ep%proj(iatom))) cycle

              do idir = 1, space%dim
                force_psi(idir) = - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
                  R_REAL(X(projector_matrix_element)(hm%ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, psi, grad_psi(:, idir, :)))
              end do

              ! We convert the forces to Cartesian coordinates
              if (space%is_periodic() .and. gr%sb%latt%nonorthogonal ) then
                force_psi(1:space%dim) = matmul(gr%sb%latt%klattice_primitive(1:space%dim, 1:space%dim), force_psi(1:space%dim))
              end if

              force_nl(1:space%dim, iatom) = force_nl(1:space%dim, iatom) + force_psi(1:space%dim)
            end do

          end if

        end do

      end if

      !The Hubbard forces
      call X(lda_u_force)(hm%lda_u, namespace, space, gr%mesh, st, iq, psib, grad_psib, force_u, allocated(hm%hm_base%phase))

      call psib%end()
      do idir = 1, space%dim
        call grad_psib(idir)%end()
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad_psi)

  ! in this case we need to convert to Cartesian coordinates at the end
  ! TODO: integrate this to the routine X(hamiltonian_elec_base_nlocal_force)
  if(hm%hm_base%apply_projector_matrices .and. .not. accel_is_enabled() .and. &
    .not. (st%symmetrize_density .and. hm%kpoints%use_symmetries)) then
    ! We convert the forces to Cartesian coordinates
    if (space%is_periodic() .and. gr%sb%latt%nonorthogonal) then
      do iatom = 1, ions%natoms
        force_nl(1:space%dim, iatom) = matmul(gr%sb%latt%klattice_primitive(1:space%dim, 1:space%dim), force_nl(1:space%dim,iatom))
      end do
    end if
  end if
 
#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, TOSTRING(X(FORCES_COMM)))
    call comm_allreduce(st%st_kpt_mpi_grp, force_nl)
    call comm_allreduce(st%st_kpt_mpi_grp, force_u)
    call comm_allreduce(st%st_kpt_mpi_grp, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif

  ! We convert the gradient of the density to cartesian coordinates before symmetrization
  ! as the two operation do not commute
  ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
  if (space%is_periodic() .and. gr%sb%latt%nonorthogonal)  then
    do ip = 1, gr%mesh%np
      grad_rho(ip, 1:space%dim) = matmul(gr%sb%latt%klattice_primitive(1:space%dim, 1:space%dim), grad_rho(ip, 1:space%dim))
    end do
  end if

  if(st%symmetrize_density) then
    call symmetrizer_init(symmetrizer, gr%mesh, gr%symm)
    SAFE_ALLOCATE(symmtmp(1:gr%mesh%np, 1:space%dim))

    call dsymmetrizer_apply(symmetrizer, gr%mesh, field_vector = grad_rho, &
              symmfield_vector = symmtmp, suppress_warning = .true.)
    grad_rho(1:gr%mesh%np, 1:space%dim) = symmtmp(1:gr%mesh%np, 1:space%dim)

    SAFE_DEALLOCATE_A(symmtmp)
    call symmetrizer_end(symmetrizer)
  end if

  call dforces_from_local_potential(gr%mesh, namespace, ions, hm%ep, grad_rho, force_loc)

  do iatom = 1, ions%natoms
    do idir = 1, space%dim
      force(idir, iatom) = force_nl(idir, iatom) + force_loc(idir, iatom) + force_u(idir, iatom)
    end do
  end do


  SAFE_DEALLOCATE_A(force_psi)
  SAFE_DEALLOCATE_A(grad_rho)

  call profiling_out(prof)
  
  POP_SUB(X(forces_from_potential))
end subroutine X(forces_from_potential)

!---------------------------------------------------------------------------
subroutine X(total_force_from_potential)(space, gr, ions, ep, st, kpoints, x, lda_u_level)
  type(space_t),                  intent(in)    :: space
  type(grid_t),                   intent(in)    :: gr
  type(ions_t),                   intent(in)    :: ions
  type(epot_t),                   intent(in)    :: ep
  type(states_elec_t),            intent(in)    :: st
  type(kpoints_t),                intent(in)    :: kpoints
  FLOAT,                          intent(inout) :: x(1:MAX_DIM)
  integer,                        intent(in)    :: lda_u_level
 
  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  FLOAT,  allocatable :: grad_rho(:, :), force(:, :)
#ifdef R_TCOMPLEX
  CMPLX :: phase
#endif

  PUSH_SUB(X(total_force_from_potential))

  ASSERT(.not. st%symmetrize_density)
  ASSERT(lda_u_level == DFT_U_NONE)

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_psi(1:np, 1:space%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:space%dim))
  grad_rho = M_ZERO
  SAFE_ALLOCATE(force(1:space%dim, 1:ions%natoms))
  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    ikpoint = st%d%get_kpoint_index(iq)
    do ist = st%st_start, st%st_end

      ff = st%d%kweights(iq) * st%occ(ist, iq) * M_TWO
      if(abs(ff) <= M_EPSILON) cycle

      call states_elec_get_state(st, gr%mesh, ist, iq, psi)

      do idim = 1, st%d%dim
        call boundaries_set(gr%der%boundaries, psi(:, idim))

        if (space%is_periodic() .and. .not. kpoints_point_is_gamma(kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:space%dim) = kpoints%get_point(ikpoint)

          !Note this phase is not correct in general. We should use the phase from the Hamiltonian
          !Here we recompute it, and moreover the vector potential is missing
#ifdef R_TCOMPLEX
          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:space%dim)*gr%mesh%x(ip, 1:space%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
          end do
#else
          ! Phase can only be applied to complex wavefunctions
          ASSERT(.false.)
#endif
        end if

        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)

        do idir = 1, space%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff*R_REAL(R_CONJ(psi(ip, idim))*grad_psi(ip, idir, idim))
          end do
        end do

      end do

      call profiling_count_operations(np*st%d%dim*space%dim*(2 + R_MUL))

      ! iterate over the projectors
      do iatom = 1, ions%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, space%dim

          force(idir, iatom) = force(idir, iatom) - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
            R_REAL(X(projector_matrix_element)(ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, psi, grad_psi(:, idir, :)))

        end do
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad_psi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, TOSTRING(X(FORCES_COMM)))
    call comm_allreduce(st%st_kpt_mpi_grp, force)
    call comm_allreduce(st%st_kpt_mpi_grp, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif

  call total_force_from_local_potential(gr%mesh, space, ep, grad_rho, x)

  do iatom = 1, ions%natoms
    do idir = 1, space%dim
      x(idir) = x(idir) - force(idir, iatom)
    end do
  end do

  SAFE_DEALLOCATE_A(force)
  POP_SUB(X(total_force_from_potential))
end subroutine X(total_force_from_potential)


! --------------------------------------------------------------------------------
subroutine X(forces_derivative)(gr, namespace, space, ions, ep, st, kpoints, lr, lr2, force_deriv, lda_u_level)
  type(grid_t),                   intent(in)    :: gr
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(ions_t),                   intent(in)    :: ions
  type(epot_t),                   intent(in)    :: ep
  type(states_elec_t),            intent(in)    :: st
  type(kpoints_t),                intent(in)    :: kpoints
  type(lr_t),                     intent(in)    :: lr
  type(lr_t),                     intent(in)    :: lr2
  CMPLX,                          intent(out)   :: force_deriv(:,:) !< (space%dim, ions%natoms)
  integer,                        intent(in)    :: lda_u_level

  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: dl_psi(:, :)
  R_TYPE, allocatable :: dl_psi2(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi2(:, :, :)
  CMPLX,  allocatable :: grad_rho(:, :)
#ifdef R_TCOMPLEX
  CMPLX :: phase
#endif
  CMPLX, allocatable  :: force_local(:, :)

  PUSH_SUB(X(forces_derivative))

  ASSERT(lda_u_level == DFT_U_NONE)

  np      = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_dl_psi(1:np, 1:space%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_dl_psi2(1:np, 1:space%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_psi(1:np, 1:space%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:space%dim))
  grad_rho = M_ZERO
  force_deriv = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi2(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    ikpoint = st%d%get_kpoint_index(iq)
    do ist = st%st_start, st%st_end

      ff = st%d%kweights(iq) * st%occ(ist, iq)
      if(abs(ff) <= M_EPSILON) cycle

      do idim = 1, st%d%dim

        call states_elec_get_state(st, gr%mesh, idim, ist, iq, psi(:, idim))
        call boundaries_set(gr%der%boundaries, psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr%X(dl_psi)(:, idim, ist, iq), dl_psi(:, idim))
        call boundaries_set(gr%der%boundaries, dl_psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr2%X(dl_psi)(:, idim, ist, iq), dl_psi2(:, idim))
        call boundaries_set(gr%der%boundaries, dl_psi2(:, idim))

        if (space%is_periodic() .and. .not. kpoints_point_is_gamma(kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:space%dim) = kpoints%get_point(ikpoint)

          !Note this phase is not correct in general. We should use the phase from the Hamiltonian
          !Here we recompute it, and moreover the vector potential is missing
#ifdef R_TCOMPLEX
          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:space%dim)*gr%mesh%x(ip, 1:space%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
            dl_psi(ip, idim) = phase*dl_psi(ip, idim)
            dl_psi2(ip, idim) = phase*dl_psi2(ip, idim)
          end do
#else
          ! Phase can only be applied to complex wavefunctions
          ASSERT(.false.)
#endif
        end if
       
        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi(:, idim), grad_dl_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi2(:, idim), grad_dl_psi2(:, :, idim), set_bc = .false.)

        !accumulate to calculate the gradient of the density
        do idir = 1, space%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff * &
              (R_CONJ(grad_psi(ip, idir, idim)) * dl_psi(ip, idim) + R_CONJ(psi(ip, idim)) * grad_dl_psi(ip, idir, idim) &
              + R_CONJ(dl_psi2(ip, idim)) * grad_psi(ip, idir, idim) + R_CONJ(grad_dl_psi2(ip, idir, idim)) * psi(ip, idim))
          end do
        end do
      end do

      ! iterate over the projectors
      do iatom = 1, ions%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, space%dim

          force_deriv(idir, iatom) = force_deriv(idir, iatom) - ff * &
            (X(projector_matrix_element)(ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, grad_psi(:, idir, :), dl_psi) &
            + X(projector_matrix_element)(ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, psi, grad_dl_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, dl_psi2, grad_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), gr%der%boundaries, st%d%dim, iq, grad_dl_psi2(:, idir, :), psi))
          
        end do
      end do

    end do
  end do
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(dl_psi)
  SAFE_DEALLOCATE_A(dl_psi2)
  SAFE_DEALLOCATE_A(grad_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi2)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, TOSTRING(X(FORCES_COMM)))
    call comm_allreduce(st%st_kpt_mpi_grp, force_deriv, dim = (/space%dim, ions%natoms/))
    call comm_allreduce(st%st_kpt_mpi_grp, grad_rho)
    call profiling_out(prof_comm)
  end if
#endif
  
  SAFE_ALLOCATE(force_local(1:space%dim, 1:ions%natoms))
  force_local = M_ZERO
  call zforces_from_local_potential(gr%mesh, namespace, ions, ep, grad_rho, force_local)
  force_deriv(:,:) = force_deriv(:,:) + force_local(:,:)
  SAFE_DEALLOCATE_A(force_local)
  SAFE_DEALLOCATE_A(grad_rho)

  POP_SUB(X(forces_derivative)) 
end subroutine X(forces_derivative)

! --------------------------------------------------------------------------------
!> lr, lr2 are wfns from electric perturbation; lr is for +omega, lr2 is for -omega.
!! for each atom, Z*(i,j) = dF(j)/dE(i)
subroutine X(forces_born_charges)(gr, namespace, space, ions, ep, st, kpoints, lr, lr2, born_charges, lda_u_level)
  type(grid_t),                   intent(in)    :: gr
  type(namespace_t),              intent(in)    :: namespace
  type(space_t),                  intent(in)    :: space
  type(ions_t),                   intent(in)    :: ions
  type(epot_t),                   intent(in)    :: ep
  type(states_elec_t),            intent(in)    :: st
  type(kpoints_t),                intent(in)    :: kpoints
  type(lr_t),                     intent(in)    :: lr(:)  !< (space%dim)
  type(lr_t),                     intent(in)    :: lr2(:) !< (space%dim)
  type(born_charges_t),           intent(inout) :: born_charges
  integer,                        intent(in)    :: lda_u_level

  integer :: iatom, idir
  CMPLX, allocatable :: force_deriv(:, :)

  PUSH_SUB(X(forces_born_charges))

  SAFE_ALLOCATE(force_deriv(1:ions%space%dim, 1:ions%natoms))

  do idir = 1, space%dim
    call X(forces_derivative)(gr, namespace, space, ions, ep, st, kpoints, lr(idir), lr2(idir), force_deriv, lda_u_level)
    do iatom = 1, ions%natoms
      born_charges%charge(:, idir, iatom) = force_deriv(:, iatom)
      born_charges%charge(idir, idir, iatom) = born_charges%charge(idir, idir, iatom) + species_zval(ions%atom(iatom)%species)
    end do
  end do

  SAFE_DEALLOCATE_A(force_deriv)

  do iatom = 1, ions%natoms
    call zsymmetrize_tensor_cart(gr%symm, born_charges%charge(:, :, iatom))
  end do

  POP_SUB(X(forces_born_charges))
end subroutine X(forces_born_charges)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
