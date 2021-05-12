!! Copyright (C) 2002-2018 M. Marques, A. Castro, A. Rubio, G. Bertsch, N. Tancogne-Dejean
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

! ---------------------------------------------------------

subroutine X(exchange_operator_single)(this, namespace, space, mesh, st_d, kpoints, ist, ik, psi, hpsi, rdmft)
  type(exchange_operator_t), intent(inout) :: this 
  type(namespace_t),         intent(in)    :: namespace
  type(space_t),             intent(in)    :: space
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  integer,                   intent(in)    :: ist
  integer,                   intent(in)    :: ik
  R_TYPE, contiguous,        intent(inout) :: psi(:, :)
  R_TYPE, contiguous,        intent(inout) :: hpsi(:, :)
  logical,                   intent(in)    :: rdmft

  type(wfs_elec_t) :: psib, hpsib

  PUSH_SUB(X(exchange_operator_single))

  call wfs_elec_init(psib, st_d%dim, ist, ist, psi, ik)
  call wfs_elec_init(hpsib, st_d%dim, ist, ist, hpsi, ik)

  call X(exchange_operator_apply)(this, namespace, space, mesh, st_d, kpoints, psib, hpsib, rdmft)

  call psib%end()
  call hpsib%end()

  POP_SUB(X(exchange_operator_single))
end subroutine X(exchange_operator_single)

! ---------------------------------------------------------
subroutine X(exchange_operator_apply)(this, namespace, space, mesh, st_d, kpoints, psib, hpsib, rdmft) 
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(space_t),             intent(in)    :: space
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  logical,                   intent(in)    :: rdmft

  PUSH_SUB(X(exchange_operator_apply))

  if(this%useACE) then
    call X(exchange_operator_apply_ACE)(this, mesh, st_d, psib, hpsib)
  else
    call X(exchange_operator_apply_standard)(this, namespace, space, mesh, st_d, kpoints, psib, hpsib, rdmft)
  end if

  POP_SUB(X(exchange_operator_apply))

end subroutine X(exchange_operator_apply)

! ---------------------------------------------------------

subroutine X(exchange_operator_apply_standard)(this, namespace, space, mesh, st_d, kpoints, psib, hpsib, rdmft)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(space_t),             intent(in)    :: space
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  logical,                   intent(in)    :: rdmft

  integer :: ibatch, jst, ip, idim, ik2, ib, ii, ist
  class(wfs_elec_t), pointer :: psi2b
  FLOAT :: exx_coef, ff
  R_TYPE, allocatable :: psi2(:, :), psi(:, :), hpsi(:, :)
  R_TYPE, allocatable :: rho(:), pot(:)
  FLOAT :: qq(1:MAX_DIM) 
  integer :: ikpoint, ikpoint2, npath
  type(fourier_space_op_t) :: coulb
  logical :: use_external_kernel

  type(profile_t), save :: prof, prof2

  PUSH_SUB(X(exchange_operator_apply_standard))

  ASSERT(associated(this%st))

  ! In case of k-points, the poisson solver must contains k-q 
  ! in the Coulomb potential, and must be changed for each q point
  exx_coef = max(this%cam_alpha,this%cam_beta)

  npath = SIZE(kpoints%coord_along_path)

  if(this%cam_beta > M_EPSILON) then
    ASSERT(this%cam_alpha < M_EPSILON)
  end if

  !The symmetries require a full treatment
  if(kpoints%use_symmetries) then
   call messages_not_implemented("symmetries with Fock operator", namespace=namespace)
  end if

  SAFE_ALLOCATE(psi(1:mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(rho(1:mesh%np))
  SAFE_ALLOCATE(pot(1:mesh%np))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st_d%dim))

  ikpoint = st_d%get_kpoint_index(psib%ik)

  use_external_kernel = (st_d%nik > st_d%spin_channels .or. this%cam_omega > M_EPSILON)
  if(use_external_kernel) then
    qq = M_ZERO
    call poisson_build_kernel(this%psolver, namespace, space, coulb, qq, this%cam_omega)
  end if


  do ibatch = 1, psib%nst
    ist = psib%ist(ibatch)
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib, ibatch, mesh%np, hpsi)

    do ik2 = 1, st_d%nik
      if(st_d%get_spin_index(ik2) /= st_d%get_spin_index(psib%ik)) cycle

      ikpoint2 = st_d%get_kpoint_index(ik2)
      !Down-sampling and q-grid
      if(st_d%nik > st_d%spin_channels) then
        if(.not.kpoints_is_compatible_downsampling(kpoints, ikpoint, ikpoint2)) cycle
        qq(1:space%dim) = kpoints%get_point(ikpoint, absolute_coordinates=.false.) &
          - kpoints%get_point(ikpoint2, absolute_coordinates=.false.)
      end if
      ! Updating of the poisson solver
      ! In case of k-points, the poisson solver must contains k-q
      ! in the Coulomb potential, and must be changed for each q point
      if(use_external_kernel) then
        call poisson_build_kernel(this%psolver, namespace, space, coulb, qq, this%cam_omega, &
          -(kpoints%full%npoints - npath)*kpoints%latt%rcell_volume*(this%singul%Fk(ik2) - this%singul%FF))
      end if

      
      do ib = 1, this%st%group%nblocks
        !We copy data into psi2b from the corresponding MPI task
        call states_elec_parallel_get_block(this%st, mesh, ib, ik2, psi2b)

        do ii = 1, psi2b%nst

          jst = psi2b%ist(ii)

          if ( .not. rdmft ) then
            ff = this%st%occ(jst, ik2)
            if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff
          else ! RDMFT
            ff = sqrt(this%st%occ(ist, psib%ik)*this%st%occ(jst, ik2)) ! Mueller functional
          end if
          ff = st_d%kweights(ik2)*exx_coef*ff

          if(ff < M_EPSILON) cycle

          call batch_get_state(psi2b, ii, mesh%np, psi2)

          call profiling_in(prof, TOSTRING(X(CODENSITIES)))
          rho = R_TOTYPE(M_ZERO)          !We compute rho_ij
          pot = R_TOTYPE(M_ZERO)

          do idim = 1, st_d%dim
            do ip = 1, mesh%np
              rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
            end do
          end do
          call profiling_out(prof)

          !and V_ij
          if(use_external_kernel) then
            call X(poisson_solve)(this%psolver, pot, rho, all_nodes = .false., kernel=coulb)
          else
            call X(poisson_solve)(this%psolver, pot, rho, all_nodes = .false.)
          end if

          !Accumulate the result
          call profiling_in(prof2, TOSTRING(X(EXCHANGE_ACCUMULATE)))
          do idim = 1, st_d%dim
            do ip = 1, mesh%np
              hpsi(ip, idim) = hpsi(ip, idim) - ff*psi2(ip, idim)*pot(ip)
            end do
          end do 
          call profiling_out(prof2)

        end do

        call states_elec_parallel_release_block(this%st, ib, psi2b)

      end do
    end do
    call batch_set_state(hpsib, ibatch, mesh%np, hpsi)

  end do

  if(use_external_kernel) then
    call fourier_space_op_end(coulb)
  end if

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator_apply_standard))
end subroutine X(exchange_operator_apply_standard)

! ---------------------------------------------------------

subroutine X(exchange_operator_apply_ACE)(this, mesh, st_d, psib, hpsib)
  type(exchange_operator_t), intent(in)    :: this
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(wfs_elec_t),          intent(in)    :: psib
  type(wfs_elec_t),          intent(inout) :: hpsib

  integer :: ibatch, ist
  FLOAT :: exx_coef
  R_TYPE, allocatable ::  psi(:, :), hpsi(:, :)
  R_TYPE :: dot

  type(profile_t), save :: prof 

  PUSH_SUB(X(exchange_operator_apply_ACE))

  exx_coef = max(this%cam_alpha,this%cam_beta)

  if(exx_coef < CNST(1.0e-3)) then
    POP_SUB(X(exchange_operator_apply_ACE))
    return
  end if

  if(.not.allocated(this%ace%X(chi))) then
    POP_SUB(X(exchange_operator_apply_ACE))
    return
  end if

  call profiling_in(prof, "EXCHANGE_APPLY_ACE")

  if(psib%is_packed() .or. hpsib%is_packed()) then
    ASSERT(psib%is_packed() .eqv. hpsib%is_packed())
  end if

  select case(psib%status())
  case(BATCH_DEVICE_PACKED, BATCH_PACKED)

    SAFE_ALLOCATE(psi(1:mesh%np, 1:st_d%dim))
    SAFE_ALLOCATE(hpsi(1:mesh%np, 1:st_d%dim))

    do ibatch = 1, psib%nst
      call batch_get_state(psib, ibatch, mesh%np, psi)
      call batch_get_state(hpsib, ibatch, mesh%np, hpsi)

      do ist = 1, this%ace%nst
        dot = X(mf_dotp)(mesh, st_d%dim, this%ace%X(chi)(:, :, ist, psib%ik), psi)
        call lalg_axpy(mesh%np, st_d%dim, -dot, this%ace%X(chi)(:, :, ist, psib%ik), hpsi)
      end do

      call batch_set_state(hpsib, ibatch, mesh%np, hpsi)
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

  case(BATCH_NOT_PACKED)

    do ibatch = 1, psib%nst
      do ist = 1, this%ace%nst
        dot = X(mf_dotp)(mesh, st_d%dim, this%ace%X(chi)(:,:,ist,psib%ik), psib%X(ff)(1:mesh%np, 1:st_d%dim, ibatch))
        call lalg_axpy(mesh%np, st_d%dim, -dot, this%ace%X(chi)(:,:, ist,psib%ik), hpsib%X(ff)(1:mesh%np, 1:st_d%dim, ibatch))
      end do
    end do

  end select


 
  call profiling_out(prof)

  POP_SUB(X(exchange_operator_apply_ACE))
end subroutine X(exchange_operator_apply_ACE)

! ---------------------------------------------------------

subroutine X(exchange_operator_compute_potentials)(this, namespace, space, mesh, st, xst, kpoints, ex, F_out)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(space_t),             intent(in)    :: space
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_t),       intent(in)    :: st 
  type(states_elec_t),       intent(inout) :: xst
  type(kpoints_t),           intent(in)    :: kpoints
  FLOAT, optional,           intent(out)   :: ex
  R_TYPE, optional,          intent(out)   :: F_out(:,:,:,:,:) !< For RDMFT

  integer :: ib, ii, ik, ist, ikloc, node_fr, node_to
  integer :: ip, idim, is, nsend, nreceiv
  integer :: kpt_start, kpt_end, st_start, st_end
  integer :: ncom, ncom_full

  R_TYPE, allocatable :: psi2(:, :), xpsi(:,:), send_buffer(:,:,:)
  R_TYPE, allocatable, target :: rec_buffer(:,:,:)
  R_TYPE, allocatable :: xpsi_ret(:,:,:), xpsi_rec(:,:,:)
  FLOAT :: exx_coef
  integer :: ikpoint, ikpoint2
  type(profile_t), save :: prof_full, prof_acc
  logical :: double_sided_communication
  type(symmetrizer_t) :: symmetrizer

#if defined(HAVE_MPI)
  integer :: send_req, status(MPI_STATUS_SIZE), mpi_err
  integer :: icom, istloc
  type(profile_t), save :: prof_comm
#endif

  PUSH_SUB(X(exchange_operator_compute_potentials))

  if(present(ex)) ex = M_ZERO

  if(present(F_out)) then
    ASSERT(.not.(st%parallel_in_states .or. st%d%kpt%parallel))
  end if
   
  !Weight of the exchange operator
  exx_coef = max(this%cam_alpha,this%cam_beta)
  if(exx_coef < CNST(1.0e-3)) then
    POP_SUB(X(exchange_operator_compute_potentials))
    return
  end if

  if(this%cam_beta > M_EPSILON) then
    ASSERT(this%cam_alpha < M_EPSILON)
  end if

  call profiling_in(prof_full, 'EXCHANGE_POTENTIALS') 

  !We deactivate the double-sided communication in case of k-points with symmetries.
  !Otherwise, we should introduce the Jacobian in doing the symmmetrization of the returned
  !This is unclear is one would gain a lot doing that
  !potential, plus some communications
  double_sided_communication = .true.
  if(kpoints%use_symmetries) double_sided_communication = .false.
  !In case of unscreened hybrids, we cannot use the double-sided communication, as 
  !the Coulomb singularity is different for k or kp
  if(this%cam_omega <= M_EPSILON) double_sided_communication = .false.

  if(kpoints%use_symmetries .and. .not. st%symmetrize_density) then
    call messages_not_implemented("ACE with KPointsUseSymmetries=yes and with SymmetrizeDensity=no")
  end if

  if(double_sided_communication) then
    !The MPI distribution scheme is not compatible with the downsampling as implemented now
    !This is because of the "returned potential"
    if(any(kpoints%downsampling(1:space%dim) /= 1)) then
      call messages_not_implemented("ACE with downsampling")
    end if
  end if

  if(kpoints%use_symmetries) then
    if(any(kpoints%downsampling(1:space%dim) /= 1)) then
      call messages_not_implemented("Downsampling with k-point symmetries")
    end if
  end if

  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(xpsi(1:mesh%np, 1:st%d%dim))

  if(.not.xst%group%block_initialized) then
    call states_elec_copy(xst, st)
  end if

  !We set to zero before doing the accumulation
  call states_elec_set_zero(xst)

  !We do the symmetrization on the received states, to reduce the amount of Poisson equation solved
  if(kpoints%use_symmetries) then
    call symmetrizer_init(symmetrizer, mesh, kpoints%symm)
  end if

  !We start by the contribution from states all present in memory
  do ik = st%d%kpt%start, st%d%kpt%end
    ikpoint = st%d%get_kpoint_index(ik)
    do ib = st%group%block_start, st%group%block_end
      !We treat a batch of states at the same time
      st_start =  st%group%block_range(ib, 1)
      st_end = st%group%block_range(ib, 2)
      SAFE_ALLOCATE(rec_buffer(1:mesh%np, 1:st%d%dim, st_start:st_end))

      do  ii = 1, st%group%psib(ib, ik)%nst
        ist =  st%group%psib(ib, ik)%ist(ii)
        call batch_get_state(st%group%psib(ib, ik), ii, mesh%np, rec_buffer(:,:,ist))
      end do

      call local_contribution( .true. )
     
      SAFE_DEALLOCATE_A(rec_buffer)
    end do !ib
  end do !ik

  !If we are not parallelized in kpt or st we are done
  if(.not.(st%parallel_in_states .or. st%d%kpt%parallel)) then
    call profiling_out(prof_full)
    SAFE_DEALLOCATE_A(psi2)
    SAFE_DEALLOCATE_A(xpsi)
    POP_SUB(X(exchange_operator_compute_potentials))
    return
  end if

  !We send to the p+s task and receive from the p-s task
  !This follows the scheme of Cardfiff et al. J. Phys.: Condens. Matter 30 095901 (2018)
  if(double_sided_communication) then
    ncom_full = int((st%st_kpt_mpi_grp%size+2)/2)-1
  else
    ncom_full = st%st_kpt_mpi_grp%size-1
  end if

  do is = 1, ncom_full
    ! node where to send the wavefunctions
    node_to = mod(st%st_kpt_mpi_grp%rank + is, st%st_kpt_mpi_grp%size)
    ! node from which we receive the wavefunctions
    node_fr = st%st_kpt_mpi_grp%rank - is
    if(node_fr < 0) node_fr = node_fr + st%st_kpt_mpi_grp%size
    node_fr = mod(node_fr, st%st_kpt_mpi_grp%size)

    if(double_sided_communication) then
      if(is == int((st%st_kpt_mpi_grp%size+2)/2)-1 &
           .and. mod(st%st_kpt_mpi_grp%size, 2) == 0 .and. node_to < st%st_kpt_mpi_grp%size/2) then
        node_to = -1
      end if

      if(is == int((st%st_kpt_mpi_grp%size+2)/2)-1 .and. mod(st%st_kpt_mpi_grp%size, 2) == 0 &
           .and. st%st_kpt_mpi_grp%rank < st%st_kpt_mpi_grp%size/2) then
        node_fr = -1
      end if
    end if

    !What we receiv for the wfn
    if(node_fr > -1) then
      st_start   = st%st_kpt_task(node_fr,1)
      st_end     = st%st_kpt_task(node_fr,2)
      kpt_start  = st%st_kpt_task(node_fr,3)
      kpt_end    = st%st_kpt_task(node_fr,4)
      nreceiv = (kpt_end-kpt_start+1) !*(st_end-st_start+1)
    else
      nreceiv = 0
    end if


    if(node_to > -1) then
      nsend = (st%d%kpt%end-st%d%kpt%start+1) !*(st%st_end-st%st_start+1)
    else
      nsend = 0
    end if

    ikloc = st%d%kpt%start
    if(node_fr > -1) then 
      ik = kpt_start
    else
      ik = -1
    end if

    !An array to store the potentials we will receive
    SAFE_ALLOCATE(send_buffer(1:mesh%np, 1:st%d%dim, st%st_start:st%st_end))
    SAFE_ALLOCATE(rec_buffer(1:mesh%np, 1:st%d%dim, st_start:st_end))
    if(double_sided_communication) then
      SAFE_ALLOCATE(xpsi_ret(1:mesh%np, 1:st%d%dim, st_start:st_end))
      SAFE_ALLOCATE(xpsi_rec(1:mesh%np, 1:st%d%dim, st%st_start:st%st_end))
    end if


    !Number of communications
    ncom = max(nsend, nreceiv)

#if defined(HAVE_MPI)
    do icom = 1, ncom
      send_req = 0

      call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')

      !Here use subarray
      ! see https://events.prace-ri.eu/event/176/contributions/57/attachments/148/296/Advanced_MPI_I.pdf

      !Sending local wfn to node_to
      if(icom<= nsend .and. node_to > -1) then
        do istloc = st%st_start, st%st_end
          call states_elec_get_state(st, mesh, istloc, ikloc, send_buffer(:,:, istloc))
        end do
        send_req = 0
        call MPI_ISend(send_buffer, mesh%np*st%d%dim*(st%st_end-st%st_start+1), R_MPITYPE, &
            node_to, icom, st%st_kpt_mpi_grp%comm, send_req, mpi_err) 
      end if

      !Receiving a wf from node_fr
      if(icom<=nreceiv .and. node_fr > -1) then
        call MPI_Recv(rec_buffer, mesh%np*st%d%dim*(st_end-st_start+1), R_MPITYPE, &
              node_fr, icom, st%st_kpt_mpi_grp%comm, status, mpi_err)
      end if

      call profiling_out(prof_comm)

      !For each wf we received, we compute the potentials 
      !using the local wfn
      if(icom<=nreceiv .and. node_fr > -1) then
        ikpoint = st%d%get_kpoint_index(ik)
        call local_contribution( .not. double_sided_communication  )
      end if

      if(send_req /= 0 ) then
        call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')
        call MPI_Wait(send_req, status, mpi_err)
        call profiling_out(prof_comm)
      end if
      send_req = 0

      !and we send it back to the task p-is (node_fr)
      if(icom<=nreceiv .and. node_fr > -1) then
        if(double_sided_communication) then
          call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')          
          send_req = 0
          call MPI_ISend(xpsi_ret, mesh%np*st%d%dim*(st_end-st_start+1), R_MPITYPE, &
            node_fr, icom, st%st_kpt_mpi_grp%comm, send_req, mpi_err)
          call profiling_out(prof_comm)
        end if

        ik = ik+1
      end if

      if(icom<=nsend .and. node_to > -1) then
        if(double_sided_communication) then
          call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')
          call MPI_Recv(xpsi_rec, mesh%np*st%d%dim*(st%st_end-st%st_start+1), R_MPITYPE, &
              node_to, icom, st%st_kpt_mpi_grp%comm, status, mpi_err)
          call profiling_out(prof_comm)

          !Accumulate the result from what was computed on the task p+s
          call profiling_in(prof_acc, "EXCHANGE_ACCUMULATE")
          !We get the potential corresponding to what we sent
          ! i.e. the wavefunction (istloc, ikloc)
          do istloc = st%st_start, st%st_end
            do idim = 1, st%d%dim 
              !TODO: Create a routine batch_add_to_state taking ist and ik as arguments
              call batch_get_state(xst%group%psib(st%group%iblock(istloc, ikloc), ikloc), &
                                        (/istloc, idim/), mesh%np, xpsi(:,idim))
              call lalg_axpy(mesh%np, M_ONE, xpsi_rec(:,idim, istloc), xpsi(:,idim))
              call batch_set_state(xst%group%psib(st%group%iblock(istloc, ikloc), ikloc), &
                                        (/istloc, idim/), mesh%np, xpsi(:,idim))
            end do
          end do

          call profiling_out(prof_acc)
        end if

        !We increment the global index
        ikloc = ikloc+1
      end if
      if(send_req /= 0 .and. double_sided_communication) then
        call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')
        call MPI_Wait(send_req, status, mpi_err)
        call profiling_out(prof_comm)
      end if
      send_req = 0

    end do
#endif
  SAFE_DEALLOCATE_A(xpsi_ret)
  SAFE_DEALLOCATE_A(xpsi_rec)
  SAFE_DEALLOCATE_A(send_buffer)
  SAFE_DEALLOCATE_A(rec_buffer)
 
  end do !is

  if(present(ex) .and. (st%parallel_in_states .or. st%d%kpt%parallel)) then
    call comm_allreduce(st%st_kpt_mpi_grp, ex)
  end if

  if(st%symmetrize_density) then
    call symmetrizer_end(symmetrizer)
  end if

  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(xpsi)

  call profiling_out(prof_full)

  POP_SUB(X(exchange_operator_compute_potentials))

  contains

  subroutine local_contribution( local )
    logical,    intent(in) :: local

    integer :: ib2, ii2, ist2, ik2
    type(profile_t), save :: prof_full
    FLOAT :: ff, ff2, kpt1(1:MAX_DIM), kpt2(1:MAX_DIM), qq(1:MAX_DIM)
    R_TYPE, allocatable :: pot(:,:), rho(:,:), ff_psi_sym(:,:)
    R_TYPE, pointer :: psi_sym(:,:), psi_sym_conj(:,:)
    integer :: iop, npath
    type(fourier_space_op_t) :: coulb
    logical :: use_external_kernel

    PUSH_SUB(X(exchange_operator_compute_potentials).local_contribution)

    call profiling_in(prof_full, "EXCHANGE_LOCAL")
     
    npath = SIZE(kpoints%coord_along_path)
    qq(1:space%dim) = M_ZERO

    use_external_kernel = (st%d%nik > st%d%spin_channels .or. this%cam_omega > M_EPSILON)
    if(use_external_kernel) then
      call poisson_build_kernel(this%psolver, namespace, space, coulb, qq, this%cam_omega)
    end if

    SAFE_ALLOCATE(pot(mesh%np, maxval(st%group%block_size)))
    SAFE_ALLOCATE(rho(mesh%np, maxval(st%group%block_size)))
    !We do not reinitialize the potential every time.
    !It is either erased completely or the previous one serves as a guess 
    pot = R_TOTYPE(M_ZERO)
    !We initialize the return buffer, if needed.
    if(.not. local) xpsi_ret = R_TOTYPE(M_ZERO)

    SAFE_ALLOCATE(ff_psi_sym(1:mesh%np, 1:st%d%dim))

    do ii = 1, kpoints_get_num_symmetry_ops(kpoints, ikpoint)
      iop = kpoints_get_symmetry_ops(kpoints, ikpoint, ii)

      if(kpoints%use_symmetries) then
        !We apply the symmetry
        kpt2(1:space%dim) = kpoints%get_point(ikpoint, absolute_coordinates=.false.)
        call symmetries_apply_kpoint_red(kpoints%symm, iop, kpt2, kpt1)
      else
        kpt1(1:space%dim) = kpoints%get_point(ikpoint, absolute_coordinates=.false.)
      end if

      !Local contribution
      do ik2 = st%d%kpt%start, st%d%kpt%end
        if(st%d%get_spin_index(ik2) /= st%d%get_spin_index(ik)) cycle
        ikpoint2 = st%d%get_kpoint_index(ik2)

        !Down-sampling and q-point
        if(use_external_kernel) then
          if(.not.kpoints_is_compatible_downsampling(kpoints, ikpoint, ikpoint2)) cycle

          kpt2(1:space%dim) = kpoints%get_point(ikpoint2, absolute_coordinates=.false.)
          qq(1:space%dim) = kpt1(1:space%dim) - kpt2(1:space%dim)
        end if
        ! Updating of the poisson solver
        ! In case of k-points, the poisson solver must contains k-q 
        ! in the Coulomb potential, and must be changed for each q point
        if(use_external_kernel) then
          call poisson_build_kernel(this%psolver, namespace, space, coulb, qq, this%cam_omega, &
                  -(kpoints%full%npoints-npath)*kpoints%latt%rcell_volume  &
                     *(this%singul%Fk(ik2)-this%singul%FF))
        end if

        !We loop over the received batch of states
        do ist = st_start, st_end

          ff = st%d%kweights(ik)*exx_coef*st%occ(ist, ik)/st%smear%el_per_state
          if(local .and. ff <= M_EPSILON) cycle

          !We do the symmetrization of the received wavefunctions
          if(kpoints%use_symmetries) then
            SAFE_ALLOCATE(psi_sym(1:mesh%np, 1:st%d%dim))
            ff = ff/kpoints_get_num_symmetry_ops(kpoints, ikpoint)
            do idim = 1, st%d%dim
              call X(symmetrizer_apply_single)(symmetrizer, mesh, iop, &
                                                   rec_buffer(:,idim,ist), psi_sym(:,idim))
            end do
          else
            psi_sym => rec_buffer(:,:,ist)
          end if

          !We precalculate some quantities
#ifdef R_TCOMPLEX
          SAFE_ALLOCATE(psi_sym_conj(1:mesh%np, 1:st%d%dim))
          do idim = 1, st%d%dim
            !$omp parallel do
            do ip = 1, mesh%np
              ff_psi_sym(ip, idim)   = ff*psi_sym(ip, idim)
              psi_sym_conj(ip, idim) = conjg(psi_sym(ip, idim))
            end do
          end do
#else
          ff_psi_sym = ff*psi_sym
          psi_sym_conj => psi_sym    
#endif

         do ib2 = st%group%block_start, st%group%block_end

           !We compute rho_ij
           !This is batchified to reuse psi_sym_conj at maximum
           if((st%group%psib(ib2, ik2)%status()) /= BATCH_DEVICE_PACKED) then
             call X(mesh_batch_codensity)(mesh, st%group%psib(ib2, ik2), psi_sym_conj, rho)
           else
             do  ii2 = 1, st%group%psib(ib2, ik2)%nst
               call batch_get_state(st%group%psib(ib2, ik2), ii2, mesh%np, psi2)
               !$omp parallel do 
               do ip = 1, mesh%np
                 rho(ip, ii2) = psi_sym_conj(ip, 1)*psi2(ip, 1)
               end do
               do idim = 2, st%d%dim
                 !$omp parallel do
                 do ip = 1, mesh%np
                   rho(ip, ii2) = rho(ip, ii2) + psi_sym_conj(ip, idim)*psi2(ip, idim)
                 end do
               end do
             end do
           end if

           !and V_ij
           !TODO: batchify poisson_solve
           if(use_external_kernel) then
             do ii2 = 1, st%group%psib(ib2, ik2)%nst
               call X(poisson_solve)(this%psolver, pot(:,ii2), rho(:,ii2), all_nodes = .false., kernel=coulb)
             end do
           else
             do ii2 = 1, st%group%psib(ib2, ik2)%nst
              call X(poisson_solve)(this%psolver, pot(:,ii2), rho(:,ii2), all_nodes = .false.)
             end do
           end if
 
           !for rdmft we need the matrix elements and no summation
           if (present(F_out)) then
             do ii2 = 1, st%group%psib(ib2, ik2)%nst
               ist2 = st%group%psib(ib2, ik2)%ist(ii2)
               F_out(:, ist, ist2, ik, ik2) = pot(:,ii2)
             end do
             cycle
           end if


           !Accumulate the result into xpsi
           call profiling_in(prof_acc, "EXCHANGE_ACCUMULATE")
           if(ff > M_EPSILON) then

             select case(xst%group%psib(ib2, ik2)%status())
             case(BATCH_DEVICE_PACKED)
               !xpsi contains the application of the Fock operator to psi2
               do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                 call batch_get_state(xst%group%psib(ib2, ik2), ii2, mesh%np, xpsi)
                 do idim = 1, st%d%dim
                   !$omp parallel do
                   do ip = 1, mesh%np
                     xpsi(ip, idim) = xpsi(ip, idim) - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
                 call batch_set_state(xst%group%psib(ib2, ik2), ii2, mesh%np, xpsi)
               end do

             case(BATCH_PACKED)
               !$omp parallel do private(ii2, idim)
               do ip = 1, mesh%np
                 do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                   do idim = 1, st%d%dim 
                     xst%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)     &
                      = xst%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)  &
                          - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
               end do            

             case(BATCH_NOT_PACKED)
               do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                 do idim = 1, st%d%dim
                   !$omp parallel do  
                   do ip = 1, mesh%np
                     xst%group%psib(ib2, ik2)%X(ff)(ip, idim, ii2)     &
                       = xst%group%psib(ib2, ik2)%X(ff)(ip, idim, ii2) &
                          - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
               end do
             end select

             !Energy, if requested
             if(present(ex)) then
               do ii2 = 1, st%group%psib(ib2, ik2)%nst
                 ist2 = st%group%psib(ib2, ik2)%ist(ii2)
                 ff2 = st%d%kweights(ik2)*st%occ(ist2, ik2)
                 ex = ex - M_HALF * ff * ff2 * R_REAL(X(mf_dotp)(mesh, rho(:,ii2), pot(:,ii2)))
               end do
             end if
           end if

           !This is what needs to be used by the sending task
           if(.not. local) then
             do  ii2 = 1, st%group%psib(ib2, ik2)%nst
               ist2 = st%group%psib(ib2, ik2)%ist(ii2)
               ff2 = st%d%kweights(ik2)*exx_coef*st%occ(ist2, ik2)/st%smear%el_per_state
               if(ff2 <= M_EPSILON) cycle

               select case(st%group%psib(ib2, ik2)%status())
               case(BATCH_DEVICE_PACKED)
                 call batch_get_state(st%group%psib(ib2, ik2), ii2, mesh%np, psi2)
                 do idim = 1, st%d%dim
                   !$omp parallel do
                   do ip = 1, mesh%np
                     xpsi_ret(ip, idim, ist) = xpsi_ret(ip, idim, ist) &
                                                 - ff2*psi2(ip, idim)*R_CONJ(pot(ip, ii2))
                   end do
                 end do
               case(BATCH_PACKED)
                 do idim = 1, st%d%dim
                   !$omp parallel do
                   do ip = 1, mesh%np
                     xpsi_ret(ip, idim, ist) = xpsi_ret(ip, idim, ist) &
                       - ff2*st%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)*R_CONJ(pot(ip, ii2))
                  end do
                end do
               case(BATCH_NOT_PACKED)
                 do idim = 1, st%d%dim
                   !$omp parallel do
                   do ip = 1, mesh%np
                     xpsi_ret(ip, idim, ist) = xpsi_ret(ip, idim, ist) &
                       - ff2*st%group%psib(ib2, ik2)%X(ff)(ip, idim, ii2)*R_CONJ(pot(ip, ii2))
                   end do
                 end do
               end select
             end do
           end if
           call profiling_out(prof_acc)
         end do !ib2

         if(kpoints%use_symmetries) then
           SAFE_DEALLOCATE_P(psi_sym)
         else
           nullify(psi_sym)
         end if
#ifdef R_TCOMPLEX
         SAFE_DEALLOCATE_P(psi_sym_conj)
#else
         nullify(psi_sym_conj)
#endif
       end do !ist
     end do !ik2
   end do !ii

   if(use_external_kernel) then
    call fourier_space_op_end(coulb)
   end if


   SAFE_DEALLOCATE_A(rho)
   SAFE_DEALLOCATE_A(pot)
   SAFE_DEALLOCATE_A(ff_psi_sym)

   call profiling_out(prof_full)

   POP_SUB(X(exchange_operator_compute_potentials).local_contribution)

  end subroutine local_contribution 

end subroutine X(exchange_operator_compute_potentials)

! ---------------------------------------------------------
! We follow the procedure defined in Lin, J. Chem. Theory Comput. 2016, 12, 2242
subroutine X(exchange_operator_ACE)(this, mesh, st, xst, phase)
  type(exchange_operator_t), intent(inout) :: this
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_t),       intent(inout) :: st
  type(states_elec_t),       intent(inout) :: xst
  CMPLX, optional,           intent(in)    :: phase(:, st%d%kpt%start:)

  integer :: nst, nst2, ierr, idim
  integer :: ik, ii, ib
  logical :: bof
  FLOAT :: exx_coef
  R_TYPE, allocatable :: MM(:,:,:), psi(:,:), xpsi(:,:)
  class(wfs_elec_t), pointer :: xpsib
  type(profile_t), save :: prof

  PUSH_SUB(X(exchange_operator_ACE))

  exx_coef = max(this%cam_alpha,this%cam_beta)

  if(exx_coef < CNST(1.0e-3)) then
    POP_SUB(X(exchange_operator_ACE))
    return
  end if

  call profiling_in(prof, 'EXCHANGE_ACE')

  !We first need to count the number of occupied states
  this%ace%nst = st%nst

  SAFE_ALLOCATE(MM(1:this%ace%nst,1:this%ace%nst, st%d%kpt%start:st%d%kpt%end))
  MM = R_TOTYPE(M_ZERO)
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim)) 
  SAFE_ALLOCATE(xpsi(1:mesh%np, 1:st%d%dim))

  if(xst%parallel_in_states) then
    call states_elec_parallel_remote_access_start(xst)
  end if

  !We loop over all states and we get the X|\psi> using MPI one sided communication
  !All the task access the memory and use it to compute their part of the M matrix,
  !then we move to the next batch of X|\psi>
  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = 1, st%group%nblocks
      !We copy data into xpsib from the corresponding MPI task
      call states_elec_parallel_get_block(xst, mesh, ib, ik, xpsib)
      do ii = 1, xpsib%nst
        nst =  st%group%psib(ib, ik)%ist(ii)
        call batch_get_state(xpsib, ii, mesh%np, xpsi)

        do nst2 = nst, st%nst
          if(nst2 < st%st_start .or. nst2 > st%st_end) cycle

          !The matrix M is hermitian as M_{ij} = <\psi_i|X|\psi_j>
          call states_elec_get_state(st, mesh, nst2, ik, psi)

          !We put a minus here as we need to perform the Cholesky factorization of -M^H
          !Eq. 13 reads for complex wfns M=MBM^H, so B = M^-T
          MM(nst2, nst, ik) = -X(mf_dotp)(mesh, st%d%dim, psi, xpsi)
        end do !nst2
      end do !ii
      call states_elec_parallel_release_block(xst, ib, xpsib)
    end do !ib
  end do !ik
  SAFE_DEALLOCATE_A(psi)

  if(xst%parallel_in_states) call states_elec_parallel_remote_access_stop(xst)

  !Reduction
  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp, MM)
  end if

  do ik = st%d%kpt%start, st%d%kpt%end
    do nst = 1, this%ace%nst
      do nst2 = nst+1, this%ace%nst
        MM(nst, nst2, ik) = R_CONJ(MM(nst2, nst, ik))
      end do
    end do
 

    !Cholesky
    bof = .false.
    ! calculate the Cholesky decomposition
    ! M =U^T*U
    call lalg_cholesky(this%ace%nst, MM(:,:,ik), bof = bof, err_code = ierr)
    if(bof) then
      write(message(1),'(a,i6,a,i6)') "The cholesky for the ACE operator failed with error code ", ierr, " for ik=", ik
      call messages_fatal(1)
    end if

    !Here we construct the ACE operator
    call lalg_invert_upper_triangular(this%ace%nst, MM(:,:,ik))  
  end do

  !The number of occupied states might have changed
  SAFE_DEALLOCATE_A(this%ace%X(chi))
  SAFE_ALLOCATE(this%ace%X(chi)(1:mesh%np, 1:st%d%dim, 1:this%ace%nst, st%d%kpt%start:st%d%kpt%end))
  this%ace%X(chi)(1:mesh%np, 1:st%d%dim, 1:this%ace%nst, st%d%kpt%start:st%d%kpt%end) = R_TOTYPE(M_ZERO)

  do ik = st%d%kpt%start, st%d%kpt%end
    do nst2 = st%st_start, st%st_end
      do idim = 1, st%d%dim
        call batch_get_state(xst%group%psib(xst%group%iblock(nst2, ik), ik), (/nst2, idim/), &
                                mesh%np, xpsi(:,idim))
        do nst = 1, this%ace%nst
          !U^-1 is upper triangular
          if(nst < nst2) cycle 
          call lalg_axpy(mesh%np, MM(nst2,nst,ik), xpsi(:, idim), this%ace%X(chi)(:, idim, nst, ik))
        end do
      end do
    end do
  end do !ik

  !Reduction
  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp,this%ace%X(chi))
  end if

#ifdef R_TCOMPLEX
  if(present(phase)) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do nst = 1, this%ace%nst
        call states_elec_set_phase(st%d, this%ace%X(chi)(:, :, nst, ik), phase(1:mesh%np, ik), mesh%np, .false.)
      end do
    end do
  end if
#endif

  SAFE_DEALLOCATE_A(MM)
  SAFE_DEALLOCATE_A(xpsi)

  call profiling_out(prof)

  POP_SUB(X(exchange_operator_ACE))

 end subroutine X(exchange_operator_ACE)

! ---------------------------------------------------------

subroutine X(exchange_operator_commute_r)(this, mesh, st_d, ik, psi, gpsi)
  type(exchange_operator_t), intent(in)    :: this
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  integer,                   intent(in)    :: ik
  R_TYPE,                    intent(in)    :: psi(:, :)
  R_TYPE,                    intent(inout) :: gpsi(:, :, :)

  integer :: idim, idir, ist
  FLOAT :: exx_coef
  R_TYPE, allocatable ::  psi_r(:, :)
  R_TYPE :: dot

  type(profile_t), save :: prof 

  PUSH_SUB(X(exchange_operator_commute_r))

  exx_coef = max(this%cam_alpha,this%cam_beta)

  if(exx_coef < CNST(1.0e-3)) then
    POP_SUB(X(exchange_operator_commute_r))
    return
  end if

  if(.not.this%useACE) &
    call messages_not_implemented("[r,V_X] without ACE")

  call profiling_in(prof, "EXCHANGE_COMMUTE_R")

  SAFE_ALLOCATE(psi_r(1:mesh%np, 1:st_d%dim))

  do ist = 1, this%ace%nst
    dot = X(mf_dotp)(mesh, st_d%dim, this%ace%X(chi)(:,:,ist,ik), psi)

    do idir = 1, mesh%sb%dim
      do idim = 1, st_d%dim
        psi_r(1:mesh%np, idim) = mesh%x(1:mesh%np,idir) * this%ace%X(chi)(1:mesh%np,idim,ist,ik)
      end do 

      call lalg_axpy(mesh%np, st_d%dim, -dot, psi_r(1:mesh%np,1:st_d%dim), gpsi(1:mesh%np, idir, 1:st_d%dim))
    end do
  end do

  do idir = 1, mesh%sb%dim
    do idim = 1, st_d%dim
      psi_r(1:mesh%np, idim) = mesh%x(1:mesh%np,idir) * psi(1:mesh%np,idim)
    end do

    do ist = 1, this%ace%nst
      dot = X(mf_dotp)(mesh, st_d%dim, this%ace%X(chi)(:,:,ist,ik), psi_r)
      call lalg_axpy(mesh%np, st_d%dim, dot, this%ace%X(chi)(:,:,ist,ik), gpsi(1:mesh%np, idir, 1:st_d%dim))
    end do
  end do
  
  SAFE_DEALLOCATE_A(psi_r)

  call profiling_out(prof)

  POP_SUB(X(exchange_operator_commute_r))
end subroutine X(exchange_operator_commute_r)

! ---------------------------------------------------------

subroutine X(exchange_operator_hartree_apply) (this, namespace, mesh, st_d, kpoints, exx_coef, psib, hpsib)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(mesh_t),              intent(in)    :: mesh
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  FLOAT,                     intent(in)    :: exx_coef
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib

  integer :: ibatch, ip, idim, ik2, ist
  FLOAT   :: ff
  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :), psi(:, :), hpsi(:, :)

  PUSH_SUB(X(exchange_operator_hartree_apply))

  if(kpoints%full%npoints > st_d%ispin) then
    call messages_not_implemented("exchange operator with k-points", namespace=namespace)
  end if

  if(this%st%parallel_in_states) then
    call messages_not_implemented("exchange operator parallel in states", namespace=namespace)
  end if

  SAFE_ALLOCATE(psi(1:mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(rho(1:mesh%np))
  SAFE_ALLOCATE(pot(1:mesh%np))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st_d%dim))

  do ibatch = 1, psib%nst
    ist = psib%ist(ibatch)
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib, ibatch, mesh%np, hpsi)
    
    do ik2 = 1, st_d%nik
      if(st_d%get_spin_index(ik2) /= st_d%get_spin_index(psib%ik)) cycle

      if(this%st%occ(ist, ik2) < M_EPSILON) cycle

      pot = R_TOTYPE(M_ZERO)
      rho = R_TOTYPE(M_ZERO)

      call states_elec_get_state(this%st, mesh, ist, ik2, psi2)

      do idim = 1, this%st%d%dim
        do ip = 1, mesh%np
          rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
        end do
      end do

      call X(poisson_solve)(this%psolver, pot, rho, all_nodes = .false.)

      ff = this%st%occ(ist, ik2)
      if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff

      do idim = 1, this%st%d%dim
        do ip = 1, mesh%np
          hpsi(ip, idim) = hpsi(ip, idim) - exx_coef*ff*psi2(ip, idim)*pot(ip)
        end do
      end do

    end do

    call batch_set_state(hpsib, ibatch, mesh%np, hpsi)

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator_hartree_apply))
end subroutine X(exchange_operator_hartree_apply)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
