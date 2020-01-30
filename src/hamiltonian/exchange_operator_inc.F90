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

subroutine X(exchange_operator_single)(this, namespace, der, st_d, ist, ik, psi, hpsi, psolver, rdmft)
  type(exchange_operator_t), intent(inout) :: this 
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  integer,                   intent(in)    :: ist
  integer,                   intent(in)    :: ik
  R_TYPE,                    intent(in)    :: psi(:, :)
  R_TYPE,                    intent(inout) :: hpsi(:, :)
  type(poisson_t),           intent(in)    :: psolver
  logical,                   intent(in)    :: rdmft

  type(wfs_elec_t) :: psib, hpsib

  PUSH_SUB(X(exchange_operator_single))

  call wfs_elec_init(psib, st_d%dim, ist, ist, psi, ik)
  call wfs_elec_init(hpsib, st_d%dim, ist, ist, hpsi, ik)

  call X(exchange_operator_apply)(this, namespace, der, st_d, psib, hpsib, psolver, rdmft)

  call psib%end()
  call hpsib%end()

  POP_SUB(X(exchange_operator_single))
end subroutine X(exchange_operator_single)

! ---------------------------------------------------------

subroutine X(exchange_operator_apply)(this, namespace, der, st_d, psib, hpsib, psolver, rdmft)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  type(poisson_t),           intent(in)    :: psolver
  logical,                   intent(in)    :: rdmft

  integer :: ibatch, jst, ip, idim, ik2, ib, ii, ist
  class(wfs_elec_t), pointer :: psi2b
  FLOAT :: exx_coef, ff
  R_TYPE, allocatable :: psi2(:, :), psi(:, :), hpsi(:, :)
  R_TYPE, allocatable :: rho(:), pot(:)
  FLOAT :: qq(1:MAX_DIM) 
  integer :: ikpoint

  type(profile_t), save :: prof, prof2

  PUSH_SUB(X(exchange_operator_apply))

  ASSERT(associated(this%st))

  ! In case of k-points, the poisson solver must contains k-q 
  ! in the Coulomb potential, and must be changed for each q point
  exx_coef = max(this%cam_alpha,this%cam_beta)

  if(this%cam_omega <= M_EPSILON) then
    if(st_d%nik > st_d%ispin) then
      call messages_not_implemented("unscreened exchange operator without k-points", namespace=namespace)
    end if
  end if

  if(der%mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)

  if(this%cam_beta > M_EPSILON) then
    ASSERT(this%cam_alpha < M_EPSILON)
  end if

  !The symmetries require a full treatment
  if(der%mesh%sb%kpoints%use_symmetries) then
   call messages_not_implemented("symmetries with Fock operator", namespace=namespace)
  end if

  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:st_d%dim))
  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:st_d%dim))

  ikpoint = states_elec_dim_get_kpoint_index(st_d, psib%ik)
  qq(1:der%dim) = M_ZERO

  do ibatch = 1, psib%nst
    ist = psib%states(ibatch)%ist
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)

    do ik2 = 1, st_d%nik
      if(states_elec_dim_get_spin_index(st_d, ik2) /= states_elec_dim_get_spin_index(st_d, psib%ik)) cycle
      
      do ib = 1, this%st%group%nblocks
        !We copy data into psi2b from the corresponding MPI task
        call states_elec_parallel_get_block(this%st, der%mesh, ib, ik2, psi2b)

        do ii = 1, psi2b%nst

          jst = psi2b%states(ii)%ist

          if ( .not. rdmft ) then
            ff = this%st%occ(jst, ik2)
            if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff
          else ! RDMFT
            ff = sqrt(this%st%occ(ist, psib%ik)*this%st%occ(jst, ik2)) ! Mueller functional
          end if
          ff = st_d%kweights(ik2)*exx_coef*ff

          if(ff < M_EPSILON) cycle

          call batch_get_state(psi2b, ii, der%mesh%np, psi2)

          call profiling_in(prof, "CODENSITIES")
          rho = R_TOTYPE(M_ZERO)          !We compute rho_ij
          pot = R_TOTYPE(M_ZERO)

          do idim = 1, st_d%dim
            do ip = 1,der%mesh%np
              rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
            end do
          end do
          call profiling_out(prof)

          !and V_ij
          call X(poisson_solve)(psolver, pot, rho, all_nodes = .false.)

          !Accumulate the result
          call profiling_in(prof2, "EXCHANGE_ACCUMULATE")
          do idim = 1, st_d%dim
            forall(ip = 1:der%mesh%np)
              hpsi(ip, idim) = hpsi(ip, idim) - ff*psi2(ip, idim)*pot(ip)
            end forall
          end do 
          call profiling_out(prof2)

        end do

        call states_elec_parallel_release_block(this%st, ib, psi2b)

      end do
    end do
    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsi)

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator_apply))
end subroutine X(exchange_operator_apply)


! ---------------------------------------------------------

subroutine X(exchange_operator_hartree_apply) (this, namespace, der, st_d, exx_coef, psib, hpsib, psolver)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  FLOAT,                     intent(in)    :: exx_coef
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  type(poisson_t),           intent(in)    :: psolver

  integer :: ibatch, ip, idim, ik2, ist
  FLOAT   :: ff
  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :), psi(:, :), hpsi(:, :)

  PUSH_SUB(X(exchange_operator))

  if(der%mesh%sb%kpoints%full%npoints > st_d%ispin) then
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
    ist = psib%states(ibatch)%ist
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)
    
    do ik2 = 1, st_d%nik
      if(states_elec_dim_get_spin_index(st_d, ik2) /= states_elec_dim_get_spin_index(st_d, psib%ik)) cycle

      if(this%st%occ(ist, ik2) < M_EPSILON) cycle

      pot = R_TOTYPE(M_ZERO)
      rho = R_TOTYPE(M_ZERO)

      call states_elec_get_state(this%st, der%mesh, ist, ik2, psi2)

      do idim = 1, this%st%d%dim
        forall(ip = 1:der%mesh%np)
          rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
        end forall
      end do

      call X(poisson_solve)(psolver, pot, rho, all_nodes = .false.)

      ff = this%st%occ(ist, ik2)
      if(st_d%ispin == UNPOLARIZED) ff = M_HALF*ff

      do idim = 1, this%st%d%dim
        forall(ip = 1:der%mesh%np)
          hpsi(ip, idim) = hpsi(ip, idim) - exx_coef*ff*psi2(ip, idim)*pot(ip)
        end forall
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
subroutine X(exchange_operator_scdm_apply) (this, namespace, scdm, der, st_d, psib, hpsib, exx_coef, hartree, psolver)
  type(exchange_operator_t), intent(in)    :: this
  type(namespace_t),         intent(in)    :: namespace
  type(scdm_t),              intent(in)    :: scdm
  type(derivatives_t),       intent(in)    :: der
  type(states_elec_dim_t),   intent(in)    :: st_d
  class(wfs_elec_t),         intent(inout) :: psib
  class(wfs_elec_t),         intent(inout) :: hpsib
  FLOAT,                     intent(in)    :: exx_coef
  logical,                   intent(in)    :: hartree
  type(poisson_t),           intent(in)    :: psolver

  integer :: ist, jst, ip, idim, ik2, ibatch
  integer :: ii, jj, kk, ll, count
  FLOAT :: ff, rr(3), dist
  R_TYPE, allocatable :: rho_l(:), pot_l(:), psil(:, :), hpsil(:, :), psi(:, :), hpsi(:, :), temp_state_global(:, :)
  type(profile_t), save :: prof_exx_scdm

  PUSH_SUB(X(exchange_operator_scdm_apply))
  
  call profiling_in(prof_exx_scdm, 'SCDM_EXX_OPERATOR')

  if(der%mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)
  
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
    ist = psib%states(ibatch)%ist
    
    call batch_get_state(psib, ibatch, der%mesh%np, psil)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsil)

    if(der%mesh%parallel_in_domains) then
#ifdef HAVE_MPI
      ! the gathering is done for the domain distribution, the states are still local to the st%mpi_grp
      call vec_allgather(der%mesh%vp, psi(:, 1), psil(:, 1))
      call vec_allgather(der%mesh%vp, hpsi(:, 1), hpsil(:, 1))
#endif
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
#ifdef HAVE_MPI
      call vec_scatter(der%mesh%vp, 0, hpsil(:, 1), hpsi(:, 1))
#endif
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
