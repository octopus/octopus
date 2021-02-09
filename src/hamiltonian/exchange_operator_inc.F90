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

subroutine X(exchange_operator_single)(this, namespace, der, st_d, kpoints, ist, ik, psi, hpsi, rdmft)
  type(exchange_operator_t), intent(inout) :: this 
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
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

  call X(exchange_operator_apply)(this, namespace, der, st_d, kpoints, psib, hpsib, rdmft)

  call psib%end()
  call hpsib%end()

  POP_SUB(X(exchange_operator_single))
end subroutine X(exchange_operator_single)

! ---------------------------------------------------------
subroutine X(exchange_operator_apply)(this, namespace, der, st_d, kpoints, psib, hpsib, rdmft) 
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  logical,                   intent(in)    :: rdmft

  PUSH_SUB(X(exchange_operator_apply))

  if(this%useACE) then
    call X(exchange_operator_apply_ACE)(this, der, st_d, psib, hpsib)
  else
    call X(exchange_operator_apply_standard)(this, namespace, der, st_d, kpoints, psib, hpsib, rdmft)
  end if

  POP_SUB(X(exchange_operator_apply))

end subroutine X(exchange_operator_apply)

! ---------------------------------------------------------

subroutine X(exchange_operator_apply_standard)(this, namespace, der, st_d, kpoints, psib, hpsib, rdmft)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
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

  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:st_d%dim))

  ikpoint = states_elec_dim_get_kpoint_index(st_d, psib%ik)

  use_external_kernel = (st_d%nik > st_d%spin_channels .or. this%cam_omega > M_EPSILON)
  if(use_external_kernel) then
    call fourier_space_op_nullify(coulb)
    qq = M_ZERO
    call poisson_build_kernel(this%psolver, namespace, der%mesh%sb, coulb, qq, this%cam_omega)
  end if


  do ibatch = 1, psib%nst
    ist = psib%ist(ibatch)
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)

    do ik2 = 1, st_d%nik
      if(states_elec_dim_get_spin_index(st_d, ik2) /= states_elec_dim_get_spin_index(st_d, psib%ik)) cycle

      ikpoint2 = states_elec_dim_get_kpoint_index(st_d, ik2)
      !Down-sampling and q-grid
      if(st_d%nik > st_d%spin_channels) then
        if(.not.kpoints_is_compatible_downsampling(kpoints, ikpoint, ikpoint2)) cycle
        qq(1:der%dim) = kpoints_get_point(kpoints, ikpoint, absolute_coordinates=.false.) &
                      - kpoints_get_point(kpoints, ikpoint2, absolute_coordinates=.false.)
      end if
      ! Updating of the poisson solver
      ! In case of k-points, the poisson solver must contains k-q
      ! in the Coulomb potential, and must be changed for each q point
      if(use_external_kernel) then
        call poisson_build_kernel(this%psolver, namespace, der%mesh%sb, coulb, qq, this%cam_omega, &
                  -(kpoints%full%npoints-npath)*der%mesh%sb%rcell_volume  &
                     *(this%singul%Fk(ik2)-this%singul%FF))
      end if

      
      do ib = 1, this%st%group%nblocks
        !We copy data into psi2b from the corresponding MPI task
        call states_elec_parallel_get_block(this%st, der%mesh, ib, ik2, psi2b)

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

          call batch_get_state(psi2b, ii, der%mesh%np, psi2)

          call profiling_in(prof, TOSTRING(X(CODENSITIES)))
          rho = R_TOTYPE(M_ZERO)          !We compute rho_ij
          pot = R_TOTYPE(M_ZERO)

          do idim = 1, st_d%dim
            do ip = 1,der%mesh%np
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
            do ip = 1, der%mesh%np
              hpsi(ip, idim) = hpsi(ip, idim) - ff*psi2(ip, idim)*pot(ip)
            end do
          end do 
          call profiling_out(prof2)

        end do

        call states_elec_parallel_release_block(this%st, ib, psi2b)

      end do
    end do
    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsi)

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

subroutine X(exchange_operator_apply_ACE)(this, der, st_d, psib, hpsib)
  type(exchange_operator_t), intent(in)    :: this
  type(derivatives_t),       intent(in)    :: der
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

  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:st_d%dim))

  do ibatch = 1, psib%nst
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)

    do ist = 1, this%ace%nst
      dot = X(mf_dotp)(der%mesh, st_d%dim, this%ace%X(chi)(:, :, ist, psib%ik), psi)
      call lalg_axpy(der%mesh%np, st_d%dim, -dot, this%ace%X(chi)(:, :, ist, psib%ik), hpsi)
    end do

    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsi)
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
 
  call profiling_out(prof)

  POP_SUB(X(exchange_operator_apply_ACE))
end subroutine X(exchange_operator_apply_ACE)

! ---------------------------------------------------------

subroutine X(exchange_operator_compute_potentials)(this, namespace, der, sb, st, kpoints)
  type(exchange_operator_t), intent(inout) :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
  type(simul_box_t),         intent(in)    :: sb
  type(states_elec_t),       intent(in)    :: st 
  type(kpoints_t),           intent(in)    :: kpoints

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

  SAFE_DEALLOCATE_A(this%ace%X(chi))
   
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
    if(any(kpoints%downsampling(1:sb%dim)/=1)) then
      call messages_not_implemented("ACE with downsampling")
    end if
  end if

  if(kpoints%use_symmetries) then
    if(any(kpoints%downsampling(1:sb%dim)/=1)) then
      call messages_not_implemented("Downsampling with k-point symmetries")
    end if
  end if

  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(xpsi(1:der%mesh%np, 1:st%d%dim))

  if(.not.this%xst%group%block_initialized) then
    call states_elec_copy(this%xst, st)
  end if

  !We set to zero before doing the accumulation
  call states_elec_set_zero(this%xst)

  !We do the symmetrization on the received states, to reduce the amount of Poisson equation solved
  if(kpoints%use_symmetries) then
    call symmetrizer_init(symmetrizer, der%mesh)
  end if

  !We start by the contribution from states all present in memory
  do ik = st%d%kpt%start, st%d%kpt%end
    ikpoint = states_elec_dim_get_kpoint_index(st%d, ik)
    do ib = st%group%block_start, st%group%block_end
      !We treat a batch of states at the same time
      st_start =  st%group%block_range(ib, 1)
      st_end = st%group%block_range(ib, 2)
      SAFE_ALLOCATE(rec_buffer(1:der%mesh%np, 1:st%d%dim, st_start:st_end))

      do  ii = 1, st%group%psib(ib, ik)%nst
        ist =  st%group%psib(ib, ik)%ist(ii)
        call batch_get_state(st%group%psib(ib, ik), ii, der%mesh%np, rec_buffer(:,:,ist))
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
    SAFE_ALLOCATE(send_buffer(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end))
    SAFE_ALLOCATE(rec_buffer(1:der%mesh%np, 1:st%d%dim, st_start:st_end))
    if(double_sided_communication) then
      SAFE_ALLOCATE(xpsi_ret(1:der%mesh%np, 1:st%d%dim, st_start:st_end))
      SAFE_ALLOCATE(xpsi_rec(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end))
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
          call states_elec_get_state(st, der%mesh, istloc, ikloc, send_buffer(:,:, istloc))
        end do
        send_req = 0
        call MPI_ISend(send_buffer, der%mesh%np*st%d%dim*(st%st_end-st%st_start+1), R_MPITYPE, &
            node_to, icom, st%st_kpt_mpi_grp%comm, send_req, mpi_err) 
      end if

      !Receiving a wf from node_fr
      if(icom<=nreceiv .and. node_fr > -1) then
        call MPI_Recv(rec_buffer, der%mesh%np*st%d%dim*(st_end-st_start+1), R_MPITYPE, &
              node_fr, icom, st%st_kpt_mpi_grp%comm, status, mpi_err)
      end if

      call profiling_out(prof_comm)

      !For each wf we received, we compute the potentials 
      !using the local wfn
      if(icom<=nreceiv .and. node_fr > -1) then
        ikpoint = states_elec_dim_get_kpoint_index(st%d, ik)
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
          call MPI_ISend(xpsi_ret, der%mesh%np*st%d%dim*(st_end-st_start+1), R_MPITYPE, &
            node_fr, icom, st%st_kpt_mpi_grp%comm, send_req, mpi_err)
          call profiling_out(prof_comm)
        end if

        ik = ik+1
      end if

      if(icom<=nsend .and. node_to > -1) then
        if(double_sided_communication) then
          call profiling_in(prof_comm, 'EXCHANGE_POTENTIALS_COMM')
          call MPI_Recv(xpsi_rec, der%mesh%np*st%d%dim*(st%st_end-st%st_start+1), R_MPITYPE, &
              node_to, icom, st%st_kpt_mpi_grp%comm, status, mpi_err)
          call profiling_out(prof_comm)

          !Accumulate the result from what was computed on the task p+s
          call profiling_in(prof_acc, "EXCHANGE_ACCUMULATE")
          !We get the potential corresponding to what we sent
          ! i.e. the wavefunction (istloc, ikloc)
          do istloc = st%st_start, st%st_end
            do idim = 1, st%d%dim 
              !TODO: Create a routine batch_add_to_state taking ist and ik as arguments
              call batch_get_state(this%xst%group%psib(st%group%iblock(istloc, ikloc), ikloc), &
                                        (/istloc, idim/), der%mesh%np, xpsi(:,idim))
              call lalg_axpy(der%mesh%np, M_ONE, xpsi_rec(:,idim, istloc), xpsi(:,idim))
              call batch_set_state(this%xst%group%psib(st%group%iblock(istloc, ikloc), ikloc), &
                                        (/istloc, idim/), der%mesh%np, xpsi(:,idim))
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
    qq(1:der%dim) = M_ZERO

    use_external_kernel = (st%d%nik > st%d%spin_channels .or. this%cam_omega > M_EPSILON)
    if(use_external_kernel) then
      call fourier_space_op_nullify(coulb)
      call poisson_build_kernel(this%psolver, namespace, sb, coulb, qq, this%cam_omega)
    end if

    SAFE_ALLOCATE(pot(der%mesh%np, maxval(st%group%block_size)))
    SAFE_ALLOCATE(rho(der%mesh%np, maxval(st%group%block_size)))
    !We do not reinitialize the potential every time.
    !It is either erased completely or the previous one serves as a guess 
    pot = R_TOTYPE(M_ZERO)
    !We initialize the return buffer, if needed.
    if(.not. local) xpsi_ret = R_TOTYPE(M_ZERO)

    SAFE_ALLOCATE(ff_psi_sym(1:der%mesh%np, 1:st%d%dim))

    do ii = 1, kpoints_get_num_symmetry_ops(kpoints, ikpoint)
      iop = kpoints_get_symmetry_ops(kpoints, ikpoint, ii)

      if(kpoints%use_symmetries) then
        !We apply the symmetry
        kpt2(1:der%dim) = kpoints_get_point(kpoints, ikpoint, absolute_coordinates=.false.)
        call symmetries_apply_kpoint_red(sb%symm, iop, kpt2, kpt1)
      else
        kpt1(1:der%dim) = kpoints_get_point(kpoints, ikpoint, absolute_coordinates=.false.)
      end if

      !Local contribution
      do ik2 = st%d%kpt%start, st%d%kpt%end
        if(states_elec_dim_get_spin_index(st%d, ik2) /= states_elec_dim_get_spin_index(st%d, ik)) cycle
        ikpoint2 = states_elec_dim_get_kpoint_index(st%d, ik2)

        !Down-sampling and q-point
        if(use_external_kernel) then
          if(.not.kpoints_is_compatible_downsampling(kpoints, ikpoint, ikpoint2)) cycle

          kpt2(1:der%dim) = kpoints_get_point(kpoints, ikpoint2, absolute_coordinates=.false.)
          qq(1:der%dim) = kpt1(1:der%dim)-kpt2(1:der%dim)
        end if
        ! Updating of the poisson solver
        ! In case of k-points, the poisson solver must contains k-q 
        ! in the Coulomb potential, and must be changed for each q point
        if(use_external_kernel) then
          call poisson_build_kernel(this%psolver, namespace, sb, coulb, qq, this%cam_omega, &
                  -(kpoints%full%npoints-npath)*sb%rcell_volume  &
                     *(this%singul%Fk(ik2)-this%singul%FF))
        end if

        !We loop over the received batch of states
        do ist = st_start, st_end

          ff = st%d%kweights(ik)*exx_coef*st%occ(ist, ik)/st%smear%el_per_state
          if(local .and. ff <= M_EPSILON) cycle

          !We do the symmetrization of the received wavefunctions
          if(kpoints%use_symmetries) then
            SAFE_ALLOCATE(psi_sym(1:der%mesh%np, 1:st%d%dim))
            ff = ff/kpoints_get_num_symmetry_ops(kpoints, ikpoint)
            do idim = 1, st%d%dim
              call X(symmetrizer_apply_single)(symmetrizer, der%mesh%np, iop, rec_buffer(:,idim,ist), psi_sym(:,idim))
            end do
          else
            psi_sym => rec_buffer(:,:,ist)
          end if

          !We precalculate some quantities
#ifdef R_TCOMPLEX
          SAFE_ALLOCATE(psi_sym_conj(1:der%mesh%np, 1:st%d%dim))
          do idim = 1, st%d%dim
            do ip = 1, der%mesh%np
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
             call X(mesh_batch_codensity)(der%mesh, st%group%psib(ib2, ik2), psi_sym_conj, rho)
           else
             do  ii2 = 1, st%group%psib(ib2, ik2)%nst
               call batch_get_state(st%group%psib(ib2, ik2), ii2, der%mesh%np, psi2)
               do ip = 1, der%mesh%np
                 rho(ip, ii2) = psi_sym_conj(ip, 1)*psi2(ip, 1)
               end do
               do idim = 2, st%d%dim
                 do ip = 1, der%mesh%np
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

           !Accumulate the result into xpsi
           call profiling_in(prof_acc, "EXCHANGE_ACCUMULATE")
           if(ff > M_EPSILON) then

             select case(this%xst%group%psib(ib2, ik2)%status())
             case(BATCH_DEVICE_PACKED)
               !xpsi contains the application of the Fock operator to psi2
               do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                 call batch_get_state(this%xst%group%psib(ib2, ik2), ii2, der%mesh%np, xpsi)
                 do idim = 1, st%d%dim
                   do ip = 1, der%mesh%np
                     xpsi(ip, idim) = xpsi(ip, idim) - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
                 call batch_set_state(this%xst%group%psib(ib2, ik2), ii2, der%mesh%np, xpsi)
               end do

             case(BATCH_PACKED)
               do ip = 1, der%mesh%np
                 do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                   do idim = 1, st%d%dim 
                     this%xst%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)     &
                      = this%xst%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)  &
                          - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
               end do            

             case(BATCH_NOT_PACKED)
               do  ii2 = 1, st%group%psib(ib2, ik2)%nst
                 do idim = 1, st%d%dim
                   do ip = 1, der%mesh%np
                     this%xst%group%psib(ib2, ik2)%X(ff)(ip, idim, ii2)     &
                       = this%xst%group%psib(ib2, ik2)%X(ff)(ip, idim, ii2) &
                          - ff_psi_sym(ip, idim)*pot(ip, ii2)
                   end do
                 end do
               end do
             end select
           end if

           !This is what needs to be used by the sending task
           if(.not. local) then
             do  ii2 = 1, st%group%psib(ib2, ik2)%nst
               ist2 = st%group%psib(ib2, ik2)%ist(ii2)
               ff2 = st%d%kweights(ik2)*exx_coef*st%occ(ist2, ik2)/st%smear%el_per_state
               if(ff2 <= M_EPSILON) cycle

               select case(st%group%psib(ib2, ik2)%status())
               case(BATCH_DEVICE_PACKED)
                 call batch_get_state(st%group%psib(ib2, ik2), ii2, der%mesh%np, psi2)
                 do idim = 1, st%d%dim
                   do ip = 1, der%mesh%np
                     xpsi_ret(ip, idim, ist) = xpsi_ret(ip, idim, ist) &
                                                 - ff2*psi2(ip, idim)*R_CONJ(pot(ip, ii2))
                   end do
                 end do
               case(BATCH_PACKED)
                 do idim = 1, st%d%dim
                   do ip = 1, der%mesh%np
                     xpsi_ret(ip, idim, ist) = xpsi_ret(ip, idim, ist) &
                       - ff2*st%group%psib(ib2, ik2)%X(ff_pack)((ii2-1)*st%d%dim+idim, ip)*R_CONJ(pot(ip, ii2))
                  end do
                end do
               case(BATCH_NOT_PACKED)
                 do idim = 1, st%d%dim
                   do ip = 1, der%mesh%np
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
subroutine X(exchange_operator_ACE)(this, der, st, phase)
  type(exchange_operator_t), intent(inout) :: this
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_t),       intent(inout) :: st
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
  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:st%d%dim)) 
  SAFE_ALLOCATE(xpsi(1:der%mesh%np, 1:st%d%dim))

  if(this%xst%parallel_in_states) then
    call states_elec_parallel_remote_access_start(this%xst)
  end if

  !We loop over all states and we get the X|\psi> using MPI one sided communication
  !All the task access the memory and use it to compute their part of the M matrix,
  !then we move to the next batch of X|\psi>
  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = 1, st%group%nblocks
      !We copy data into xpsib from the corresponding MPI task
      call states_elec_parallel_get_block(this%xst, der%mesh, ib, ik, xpsib)
      do ii = 1, xpsib%nst
        nst =  st%group%psib(ib, ik)%ist(ii)
        call batch_get_state(xpsib, ii, der%mesh%np, xpsi)

        do nst2 = nst, st%nst
          if(nst2 < st%st_start .or. nst2 > st%st_end) cycle

          !The matrix M is hermitian as M_{ij} = <\psi_i|X|\psi_j>
          call states_elec_get_state(st, der%mesh, nst2, ik, psi)

          !We put a minus here as we need to perform the Cholesky factorization of -M^H
          !Eq. 13 reads for complex wfns M=MBM^H, so B = M^-T
          MM(nst2, nst, ik) = -X(mf_dotp)(der%mesh, st%d%dim, psi, xpsi)
        end do !nst2
      end do !ii
      call states_elec_parallel_release_block(this%xst, ib, xpsib)
    end do !ib
  end do !ik
  SAFE_DEALLOCATE_A(psi)

  if(this%xst%parallel_in_states) call states_elec_parallel_remote_access_stop(this%xst)

  !Reduction
  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp%comm, MM)
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
  SAFE_ALLOCATE(this%ace%X(chi)(1:der%mesh%np, 1:st%d%dim, 1:this%ace%nst, st%d%kpt%start:st%d%kpt%end))
  this%ace%X(chi)(1:der%mesh%np, 1:st%d%dim, 1:this%ace%nst, st%d%kpt%start:st%d%kpt%end) = R_TOTYPE(M_ZERO)

  do ik = st%d%kpt%start, st%d%kpt%end
    do nst2 = st%st_start, st%st_end
      do idim = 1, st%d%dim
        call batch_get_state(this%xst%group%psib(this%xst%group%iblock(nst2, ik), ik), (/nst2, idim/), &
                                der%mesh%np, xpsi(:,idim))
        do nst = 1, this%ace%nst
          !U^-1 is upper triangular
          if(nst < nst2) cycle 
          call lalg_axpy(der%mesh%np, MM(nst2,nst,ik), xpsi(:, idim), this%ace%X(chi)(:, idim, nst, ik))
        end do
      end do
    end do
  end do !ik

  !Reduction
  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp%comm,this%ace%X(chi))
  end if

#ifdef R_TCOMPLEX
  if(present(phase)) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do nst = 1, this%ace%nst
        call states_elec_set_phase(st%d, this%ace%X(chi)(:, :, nst, ik), phase(1:der%mesh%np, ik), der%mesh%np, .false.)
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

subroutine X(exchange_operator_hartree_apply) (this, namespace, der, st_d, kpoints, exx_coef, psib, hpsib)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
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

  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:st_d%dim))

  do ibatch = 1, psib%nst
    ist = psib%ist(ibatch)
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)
    
    do ik2 = 1, st_d%nik
      if(states_elec_dim_get_spin_index(st_d, ik2) /= states_elec_dim_get_spin_index(st_d, psib%ik)) cycle

      if(this%st%occ(ist, ik2) < M_EPSILON) cycle

      pot = R_TOTYPE(M_ZERO)
      rho = R_TOTYPE(M_ZERO)

      call states_elec_get_state(this%st, der%mesh, ist, ik2, psi2)

      do idim = 1, this%st%d%dim
        do ip = 1, der%mesh%np
          rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
        end do
      end do

      call X(poisson_solve)(this%psolver, pot, rho, all_nodes = .false.)

      ff = this%st%occ(ist, ik2)
      if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff

      do idim = 1, this%st%d%dim
        do ip = 1, der%mesh%np
          hpsi(ip, idim) = hpsi(ip, idim) - exx_coef*ff*psi2(ip, idim)*pot(ip)
        end do
      end do

    end do

    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsi)

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator_hartree_apply))
end subroutine X(exchange_operator_hartree_apply)

! scdm_EXX
! ---------------------------------------------------------
subroutine X(exchange_operator_scdm_apply) (this, namespace, scdm, der, st_d, kpoints, &
                                               psib, hpsib, exx_coef, hartree)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(scdm_t),              intent(in)    :: scdm
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  type(kpoints_t),           intent(in)    :: kpoints
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  FLOAT,                     intent(in)    :: exx_coef
  logical,                   intent(in)    :: hartree

  integer :: ist, jst, ip, idim, ik2, ibatch
  integer :: ii, jj, kk, ll, count
  FLOAT :: ff, rr(3), dist
  R_TYPE, allocatable :: rho_l(:), pot_l(:), psil(:, :), hpsil(:, :), psi(:, :), hpsi(:, :), temp_state_global(:, :)
  type(profile_t), save :: prof_exx_scdm

  PUSH_SUB(X(exchange_operator_scdm_apply))
  
  call profiling_in(prof_exx_scdm, TOSTRING(X(SCDM_EXX_OPERATOR)))

  if(kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)
  
  ! make sure scdm is localized
  call X(scdm_localize)(scdm, namespace, this%st, der%mesh)
  
  SAFE_ALLOCATE(psil(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsil(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(psi(1:der%mesh%np_global, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np_global, 1:st_d%dim))
  SAFE_ALLOCATE(temp_state_global(der%mesh%np_global, this%st%d%dim))
  SAFE_ALLOCATE(rho_l(1:this%scdm%full_box))
  SAFE_ALLOCATE(pot_l(1:this%scdm%full_box))
  
  do ibatch = 1, psib%nst
    ist = psib%ist(ibatch)
    
    call batch_get_state(psib, ibatch, der%mesh%np, psil)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsil)

    if(der%mesh%parallel_in_domains) then
      ! the gathering is done for the domain distribution, the states are still local to the st%mpi_grp
      call vec_allgather(der%mesh%vp, psi(:, 1), psil(:, 1))
      call vec_allgather(der%mesh%vp, hpsi(:, 1), hpsil(:, 1))
    else
      psi(1:der%mesh%np, 1:st_d%dim) = psil(1:der%mesh%np, 1:st_d%dim)
      hpsi(1:der%mesh%np, 1:st_d%dim) = hpsil(1:der%mesh%np, 1:st_d%dim)
    end if
    
    ! accumulate exchange contribution to Hpsi in a temp array and add to Hpsi at the end
    temp_state_global(:,:) = M_ZERO

    do ik2 = 1, st_d%nik
      if(states_elec_dim_get_spin_index(st_d, ik2) /= states_elec_dim_get_spin_index(st_d, psib%ik)) cycle
      count = 0
      do jst = this%scdm%st_exx_start, this%scdm%st_exx_end

        if(this%st%occ(jst, ik2) < M_EPSILON) cycle
        ! for psi in scdm representation check if it overlaps with the box of jst
        ! NOTE: this can be faster by building an array with overlapping index pairs
        !       within the radius of scdm%box_size
        if(this%scdm%psi_scdm) then
          do ii = 1, 3
            rr(1:3) = this%scdm%center(ii,ist) - this%scdm%center(ii,jst)
          end do
          dist = sqrt(dot_product(rr, rr))
          if(dist .gt. this%scdm%box_size) cycle
        end if

        ! in Hartree we just remove the self-interaction
        if(hartree .and. jst /= ist) cycle

        ! for scdm do product only in the local box
        rho_l(:) = M_ZERO

        ! copy density to local box
        do jj = 1, this%scdm%box_size*2 + 1
          do kk = 1, this%scdm%box_size*2 + 1
            do ll = 1, this%scdm%box_size*2 + 1
              ip = (jj - 1)*((this%scdm%box_size*2 + 1))**2+(kk - 1)*((this%scdm%box_size*2 + 1)) + ll
              rho_l(ip) = R_CONJ(this%scdm%X(psi)(ip, jst))*psi(this%scdm%box(jj, kk, ll, jst), 1)
            end do
          end do
        end do

        call X(poisson_solve)(this%scdm%poisson, pot_l, rho_l, all_nodes=.false.)

        ff = this%st%occ(jst, ik2)
        if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff

        do idim = 1, this%st%d%dim
          ! potential in local box to full H*psi 
          do jj =1, this%scdm%box_size*2 + 1
            do kk =1, this%scdm%box_size*2 + 1
              do ll =1, this%scdm%box_size*2 + 1
                ip = (jj - 1)*((this%scdm%box_size*2 + 1))**2 + (kk - 1)*((this%scdm%box_size*2 + 1)) + ll
                temp_state_global(this%scdm%box(jj, kk, ll, jst), idim) = &
                  temp_state_global(this%scdm%box(jj, kk, ll, jst), idim) - exx_coef*ff*this%scdm%X(psi)(ip, jst)*pot_l(ip)
              end do
            end do
          end do

        end do

      end do
    end do

    ! sum contributions to hpsi from all processes in the st_exx_grp group
    call comm_allreduce(this%scdm%st_exx_grp%comm, temp_state_global)
    
    ! add exchange contribution to the input state
    hpsi(1:der%mesh%np_global, 1) =  hpsi(1:der%mesh%np_global, 1) + temp_state_global(1:der%mesh%np_global, 1)

    if(der%mesh%parallel_in_domains) then
      call vec_scatter(der%mesh%vp, 0, hpsil(:, 1), hpsi(:, 1))
    else
      hpsil(1:der%mesh%np, 1:st_d%dim) = hpsi(1:der%mesh%np, 1:st_d%dim)
    end if

    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsil)
  end do
  
  SAFE_DEALLOCATE_A(psil)
  SAFE_DEALLOCATE_A(hpsil)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(temp_state_global)
  SAFE_DEALLOCATE_A(rho_l)
  SAFE_DEALLOCATE_A(pot_l)

  call profiling_out(prof_exx_scdm)
  
  POP_SUB(X(exchange_operator_scdm_apply))
end subroutine X(exchange_operator_scdm_apply)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
