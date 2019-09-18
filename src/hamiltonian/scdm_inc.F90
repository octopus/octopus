!! Copyright (C) 2015 H. Huebener
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

!> this performs the SCDM localization, transforming the original set of states KS (Kohn-Sham)
!! into the set SCDM, by first performing RRQR and then Cholesky for orthogonalization

subroutine X(scdm_localize)(st,mesh,scdm)

  type(states_elec_t), intent(in) :: st !< this contains the non-localize set KS (for now from hm%hf_st which is confusing)
  type(mesh_t), intent(in)   :: mesh
  type(scdm_t) :: scdm

  integer :: ii, jj, kk, ll, vv, count, ip, nval, info,i1, i2, i3, idim, j1, j2, j3
  integer, allocatable :: JPVT(:)
  integer :: icenter(3), ind_center

  integer :: nn(3)
  R_TYPE, allocatable :: SCDM_temp(:,:), Pcc(:,:), SCDM_matrix(:,:)
  FLOAT  :: error
  R_TYPE, allocatable ::  state_global(:), temp_state(:,:),  temp_column(:) !< work arrays
  integer,  allocatable :: temp_box(:,:,:)
  integer :: ix(3), nr(3,2)
  logical :: out_of_index_range(3), out_of_mesh(3)
  FLOAT, allocatable :: lxyz_domains(:,:),  lxyz_global(:), lxyz_local(:)

  PUSH_SUB(X(scdm_localize))
  ! check if already localized
  if(scdm_is_local) then
    POP_SUB(X(scdm_localize))
    return
  end if

  call profiling_in(prof_scdm,"SCDM")

  nval = st%nst ! TODO: check that this is really the number of valence states

  if (st%d%nik /= 1 .or. st%d%dim /= 1) call messages_not_implemented("SCDM with k-points or dims")

  SAFE_ALLOCATE(JPVT(1:mesh%np_global))
  call X(scdm_rrqr)(st,scdm, mesh, nval,scdm%root,1, jpvt)

  !form SCDM_matrix rows for states that are local in the scdm%st_grp
  SAFE_ALLOCATE(SCDM_matrix(1:scdm%st%lnst,nval))
  SAFE_ALLOCATE(temp_state(1:mesh%np,1))
  SAFE_ALLOCATE(state_global(1:mesh%np_global))
  count = 0
  do vv = scdm%st%st_start, scdm%st%st_end
    count = count +1
    call states_elec_get_state(st, mesh, vv, st%d%nik, temp_state)
#ifdef HAVE_MPI
    call vec_allgather(mesh%vp, state_global, temp_state(1:mesh%np,1))
#else
    state_global(1:mesh%np)=temp_state(1:mesh%np,1)
#endif
    ! loop over JPVT to copy columns of the density matrix
    do ii=1,nval
      SCDM_matrix(count,ii) = state_global(JPVT(ii))
    end do
  end do

  call profiling_in(prof_scdm_matmul1,"SCDM_matmul1")

  SAFE_ALLOCATE(SCDM_temp(1:mesh%np,1:nval))
  SAFE_ALLOCATE(temp_column(1:mesh%np))
  SCDM_temp(:,:) = M_ZERO
  do ii = 1, nval
    ! form the SCDM states in by performing the sum over
    ! SCDM elements respecting domain and state distribution
    temp_column(1:mesh%np) = M_ZERO
    count = 0
    do vv = scdm%st%st_start,scdm%st%st_end
      count = count +1
      call states_elec_get_state(st, mesh, vv, st%d%nik, temp_state)
      temp_column(1:mesh%np) = temp_column(1:mesh%np) + &
                                  temp_state(1:mesh%np,1)* R_CONJ(SCDM_matrix(count,ii))
    end do
    SCDM_temp(1:mesh%np,ii) = temp_column(1:mesh%np)
  end do
  ! the above is a prtial sum in states distribution
  call comm_allreduce(scdm%st_grp%comm, SCDM_temp, (/mesh%np, nval/))
  call profiling_out(prof_scdm_matmul1)

  SAFE_DEALLOCATE_A(temp_state)
  SAFE_DEALLOCATE_A(state_global)
  SAFE_DEALLOCATE_A(temp_column)

  ! --- Orthogoalization ----
  ! form lower triangle of Pcc
  SAFE_ALLOCATE(Pcc(1:nval,1:nval))
  Pcc(:,:) = M_ZERO
  do ii = 1, nval
    do jj = 1, ii
      do vv = 1,scdm%st%lnst
        Pcc(ii,jj) = Pcc(ii,jj)+ SCDM_matrix(vv,ii)*R_CONJ(SCDM_matrix(vv,jj))
      end do
    end do
  end do
  call comm_allreduce(scdm%st_grp%comm,Pcc)
  SAFE_DEALLOCATE_A(SCDM_matrix)

  ! Cholesky fact.
  call X(POTRF)("L", nval, Pcc, nval, info )
  if (info /= 0) then
    if (info < 0) then
      write(message(1),'(A28,I2)') 'Illegal argument in DPOTRF: ', info
      call messages_fatal(1, namespace=st%namespace)
    else
      message(1) = 'Fail of Cholesky, not pos-semi-def '
      call messages_fatal(1, namespace=st%namespace)
    end if
    stop
  end if

  ! transpose
  Pcc(:,:) = transpose(R_CONJ(Pcc(:,:)))
  ! invert
  call X(invert)(nval,Pcc)

  call profiling_in(prof_scdm_matmul3,"SCDM_matmul3")
  ! form orthogonal SCDM
  SAFE_ALLOCATE(temp_state(1:mesh%np,1))
  do vv=scdm%st%st_start,scdm%st%st_end
    temp_state(1:mesh%np,:)= M_ZERO
    do jj = 1, nval
      temp_state(1:mesh%np,1) = temp_state(1:mesh%np,1) + SCDM_temp(1:mesh%np,jj)*Pcc(jj, vv)
    end do
    call states_elec_set_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
  end do

  SAFE_DEALLOCATE_A(Pcc)
  SAFE_DEALLOCATE_A(SCDM_temp)

  call profiling_out(prof_scdm_matmul3)

  ! normalise SCDM states
  do vv = scdm%st%st_start,scdm%st%st_end
    call states_elec_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
    call X(mf_normalize)(mesh, scdm%st%d%dim, temp_state)
    call states_elec_set_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
  end do

  ! end of SCDM procedure

  ! find centers, by computing center of mass of |psi|^2
  scdm%center(:,:) = M_ZERO
  count = 0
  ! get domains of mesh%idx%lxyz
  SAFE_ALLOCATE(lxyz_domains(1:mesh%np,3))
  SAFE_ALLOCATE(lxyz_global(1:mesh%np_global))
  SAFE_ALLOCATE(lxyz_local(1:mesh%np))
  do ii=1,3
    lxyz_global(1:mesh%np_global) = mesh%idx%lxyz(1:mesh%np_global,ii)
#ifdef HAVE_MPI
    call vec_scatter(mesh%vp, 0, lxyz_local, lxyz_global)
#endif
    lxyz_domains(1:mesh%np,ii) = lxyz_local(1:mesh%np)
  end do

  do vv = scdm%st%st_start,scdm%st%st_end
    call states_elec_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
    do ii=1,3
      scdm%center(ii,vv) = sum(temp_state(1:mesh%np,1)*R_CONJ(temp_state(1:mesh%np,1))*&
           lxyz_domains(1:mesh%np,ii)*mesh%spacing(ii))*mesh%volume_element
    end do
  end do
    ! reduce to world, since the above can be state+domain distribution
  call comm_allreduce(mpi_world%comm, scdm%center)

  SAFE_DEALLOCATE_A(lxyz_domains)
  SAFE_DEALLOCATE_A(lxyz_global)
  SAFE_DEALLOCATE_A(lxyz_local)

  ! -------------- end of SCDM procedure -----------------

  ! ------------ copy local box of state ---------------------

  error = M_ZERO
  scdm%X(psi)(:,:) =  M_ZERO
  scdm%box(:,:,:,:) =  M_ZERO
  SAFE_ALLOCATE(temp_box(1:2*scdm%box_size+1,1:2*scdm%box_size+1,1:2*scdm%box_size+1))
  SAFE_ALLOCATE(state_global(1:mesh%np_global))
  do vv = scdm%st%st_start,scdm%st%st_end

    ! find integer index of center
    do ii = 1, 3
      icenter(ii) = scdm%center(ii,vv)/mesh%spacing(ii)
    end do
    ! find index of center in the mesh
    ind_center = mesh%idx%lxyz_inv(icenter(1),icenter(2),icenter(3))

    ! make sure that box does not fall out of range of the index structure
    call check_box_in_index(mesh%idx,icenter(:),scdm%box_size,out_of_index_range)

    ! only periodic dimensions can be out of range
    do idim=1,3
      if(out_of_index_range(idim).and.idim > mesh%sb%periodic_dim) then
        message(1) = 'SCDM box out of index range in non-periodic dimension'
        call messages_fatal(1, namespace=st%namespace)
      end if
    end do

    ! make list with points in the box
    if (all(out_of_index_range .eqv. (/.false.,.false.,.false./)) ) then
      temp_box(:,:,:) =  mesh%idx%lxyz_inv(icenter(1)-scdm%box_size:icenter(1)+scdm%box_size, &
           icenter(2)-scdm%box_size:icenter(2)+scdm%box_size, &
           icenter(3)-scdm%box_size:icenter(3)+scdm%box_size)

      ! check if all indices are within the mesh
      out_of_mesh(1:3) = .false.
      do idim=1,3
        if(minval(minval(temp_box,dim=idim)) < 1 .or. &
           maxval(maxval(temp_box,dim=idim)) > mesh%np_global) then
          out_of_mesh(idim) = .true.
          ! can only be out of mesh in periodic direction
          if(idim > mesh%sb%periodic_dim ) then
            message(1) = 'SCDM box out of mesh in non-periodic dimension'
            call messages_fatal(1, namespace=st%namespace)
          end if
        end if
      end do

    end if

    ! in case there are periodic replica go through every point
    ! NOTE: in principle this would be needed only for the periodic
    !       dimensions, because above we have made sure that in
    !       all other directions there are no points out of the mesh
    if(any(out_of_mesh        .eqv. (/.true.,.true.,.true./)) .or. &
       any(out_of_index_range .eqv. (/.true.,.true.,.true./)) ) then

      ! dimension of full simulation box
      nn = scdm%full_cube_n
      ! limits of the indices that are on the mesh
      ! NOTE: this should already be defined somewhere, but I didnt find it
      do idim=1,3
        if(mod(scdm%full_cube_n(idim),2) == 0 ) then
          nr(idim,1) = -nn(idim)/2
          nr(idim,2) =  nn(idim)/2-1
        else
          nr(idim,1) = -nn(idim)/2
          nr(idim,2) =  nn(idim)/2
        end if
      end do

      do  i1 = -scdm%box_size, scdm%box_size
        do i2 = -scdm%box_size, scdm%box_size
          do i3 = -scdm%box_size, scdm%box_size

            ix(:) = icenter(:)+(/i1,i2,i3/)

            do idim=1,3
              if( ix(idim) < nr(idim,1)) then
                ix(idim) = ix(idim) + nn(idim)
              else if( ix(idim) > nr(idim,2)) then
                ix(idim) = ix(idim) - nn(idim)
              end if
            end do

            ! indices of box are 1-based
            j1 = i1 + scdm%box_size + 1
            j2 = i2 + scdm%box_size + 1
            j3 = i3 + scdm%box_size + 1
            temp_box(j1,j2,j3) = mesh%idx%lxyz_inv(ix(1),ix(2),ix(3))

            !if(temp_box(j1,j2,j3) < 1 .or. temp_box(j1,j2,j3) > mesh%np_global) then
            !  print *, 'fail'
            !  print *, nr
            ! print *, ix
            !end if

          end do!i3
        end do!i2
      end do!i1

    end if

    ! check that box is well defined now
    if(minval(temp_box) <= 0.or.maxval(temp_box) > mesh%np_global ) then
      message(1) = 'SCDM box mapping failed'
      call messages_fatal(1, namespace=st%namespace)
    end if

    scdm%box(:,:,:,vv) = temp_box(:,:,:)

    ! to copy the scdm state in the box, we need the global state (this is still in state distribution)
    call states_elec_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
#ifdef HAVE_MPI
    state_global(1:mesh%np_global) = M_ZERO
    call vec_allgather(mesh%vp, state_global, temp_state(1:mesh%np,1))
#endif
    ! copy points to box
    ! this box refers to the global mesh
    do jj = 1, scdm%box_size*2+1
      do kk = 1, scdm%box_size*2+1
        do ll = 1, scdm%box_size*2+1
          ! map into the box
          ip = (jj-1)*((scdm%box_size*2+1))**2+(kk-1)*((scdm%box_size*2+1)) + ll
          scdm%X(psi)(ip,vv) = state_global(scdm%box(jj,kk,ll,vv))
        end do
      end do
    end do

    ! compute localization error
    error = error + M_ONE - dot_product(scdm%X(psi)(:,vv),scdm%X(psi)(:,vv))*mesh%volume_element

  end do

  SAFE_DEALLOCATE_A(temp_box)
  SAFE_DEALLOCATE_A(temp_state)
  SAFE_DEALLOCATE_A(state_global)

  ! the boxed SCDM states as well as their mapping boxes need to be available globally
  ! NOTE: strictly speaking only within their respective scdm%st_exx_grp, but that requires
  !       some more complicated communication, so for now this is global
  call comm_allreduce(scdm%st_grp%comm, scdm%X(psi))
  call comm_allreduce(scdm%st_grp%comm, scdm%box)

  ! sum of the error of the norm of SCDM states by truncating them
  call comm_allreduce(scdm%st_grp%comm, error)

  if (scdm%root .and. scdm%verbose) call messages_print_var_value(stdout, 'SCDM localization error:', error/st%nst)

  ! set flag to do this only once
  scdm_is_local = .true.

  call profiling_out(prof_scdm)

  POP_SUB(X(scdm_localize))

end subroutine X(scdm_localize)

!> rotate states from KS to SCDM representation and construct SCDM states
subroutine X(scdm_rotate_states)(st,mesh,scdm)
  type(states_elec_t), intent(inout)  :: st
  type(mesh_t),      intent(in)       :: mesh
  type(scdm_t),      intent(inout)    :: scdm

  PUSH_SUB(X(scdm_rotate_states))

  ! create localized SCDM representation of the states in st
  scdm_is_local = .false.
  call X(scdm_localize)(st,mesh,scdm)

  ! overwrite state object with the scdm states
  call states_elec_copy(st,scdm%st)

  POP_SUB(X(scdm_rotate_states))

end subroutine X(scdm_rotate_states)


!> stupid routine to invert with LAPACK
!! is very redundant here, shoudl be replaced by something smart
subroutine X(invert)(nn, A)
  integer         :: nn
  R_TYPE          :: A(nn,nn)

  integer         :: ierror,ipiv(nn), lwork
  R_TYPE,pointer  :: work(:)
  FLOAT           :: temp

  PUSH_SUB(X(invert))

  call X(getrf)(nn, nn, A, nn, ipiv, ierror )

  if( ierror == 0 ) then
    !workspace query
    call X(getri)(nn, A, nn, ipiv, temp, -1, ierror )
    lwork = temp ! dimension of workspace
    allocate(work(lwork*2))
    call X(getri)(nn, A, nn, ipiv, work, lwork, ierror )
  else
    message(1) = 'Terminating due to failed LU decomp'
    call messages_fatal(1)
  end if
  if (ierror /= 0) then
    message(1) = 'Terminating due to failed inversion'
    call messages_fatal(1)
  end if
  deallocate(work)

  POP_SUB(X(invert))

end subroutine X(invert)

!> Perform RRQR on the transpose states stored in the states object
!! and return the pivot vector 
!! This is not an all-purose routien for RRQR, but only operates on the
!! specific set stored in st
subroutine X(scdm_rrqr)(st, scdm, mesh, nst,root, ik, jpvt)
  type(states_elec_t), intent(in) :: st
  type(scdm_t), intent(in)        :: scdm  !< this is only needed for the proc_grid
  type(mesh_t), intent(in)        :: mesh
  integer, intent(in)             :: nst
  logical, intent(in)             :: root !< this is needed for serial 
  integer, intent(in)             :: ik ! perform SCDM with this k-point
  integer, intent(out)            :: jpvt(:)

  integer :: total_np, nref, info, wsize
  R_TYPE, allocatable :: tau(:), work(:)
  R_TYPE :: tmp
  FLOAT, allocatable :: rwork(:)
  R_TYPE, allocatable ::  state_global(:), temp_state(:,:)
  R_TYPE, allocatable :: KSt(:,:)
  R_TYPE, allocatable :: psi(:, :)
  integer :: ii,ist,  count, lnst
  logical :: do_serial
  integer :: psi_block(2), blacs_info
  integer, allocatable :: ipiv(:)
#ifdef HAVE_SCALAPACK
  integer :: psi_desc(BLACS_DLEN)
#ifndef R_TREAL
  integer :: rwsize
  FLOAT :: tmp2
#endif
#endif
#ifdef HAVE_MPI
  integer :: sender
#endif

  PUSH_SUB(X(scdm_rrqr))
  call profiling_in(prof_scdm_QR,"SCDM_QR")

  ASSERT(.not. mesh%use_curvilinear)
  ASSERT(nst == st%nst)

  lnst = st%lnst

  ! decide whether we can use ScaLAPACK
  do_serial = .false.
  if(mesh%parallel_in_domains .or. st%parallel_in_states) then
#ifndef HAVE_SCALAPACK
     message(1) = 'The RRQR is performed in serial. Try linking ScaLAPCK'
     call messages_warning(1, namespace=st%namespace)
     do_serial = .true.
#else
     if(.not.st%scalapack_compatible) then
        message(1) = 'The RRQR is performed in serial. Try setting ScaLAPACKCompatible = yes'
        call messages_warning(1, namespace=st%namespace)
        do_serial = .true.
     end if
#endif
  else
     do_serial = .true.
  endif

  if(.not.do_serial) then
    
    call states_elec_parallel_blacs_blocksize(st, mesh, psi_block, total_np)
    
    ! allocate local part of transpose state matrix
    SAFE_ALLOCATE(KSt(1:lnst,1:total_np))
    SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
    
    ! copy states into the transpose matrix
    count = 0
    do ist = st%st_start,st%st_end
      count = count + 1

      call states_elec_get_state(st, mesh, ist, ik, psi)
      
      ! We need to set to zero some extra parts of the array
      if(st%d%dim == 1) then
        psi(mesh%np + 1:psi_block(1), 1:st%d%dim) = M_ZERO
      else
        psi(mesh%np + 1:mesh%np_part, 1:st%d%dim) = M_ZERO
      end if

      KSt(count, 1:total_np) = psi(1:total_np, 1)
    end do

    SAFE_DEALLOCATE_A(psi)
     
    ! DISTRIBUTE THE MATRIX ON THE PROCESS GRID
    ! Initialize the descriptor array for the main matrices (ScaLAPACK)
#ifdef HAVE_SCALAPACK
    call descinit(psi_desc(1), nst, total_np, psi_block(2), psi_block(1), 0, 0, &
      scdm%proc_grid%context, lnst, blacs_info)
#endif
     
    if(blacs_info /= 0) then
       write(message(1),'(a,i6)') 'descinit failed with error code: ', blacs_info
       call messages_fatal(1, namespace=st%namespace)
    end if
    
    nref = min(nst, total_np)
    SAFE_ALLOCATE(tau(1:nref))
    tau = M_ZERO

    ! calculate the QR decomposition
    SAFE_ALLOCATE(ipiv(1:total_np))
    ipiv(1:total_np) = 0

    ! Note: lapack routine has different number of arguments depending on type
#ifdef HAVE_SCALAPACK
#ifndef R_TREAL
    call pzgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), tmp, -1, tmp2, -1, blacs_info) 
#else 
    call pdgeqpf( nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), tmp, -1, blacs_info)
#endif
#endif
    
    if(blacs_info /= 0) then
      write(message(1),'(a,i6)') 'scalapack geqrf workspace query failed with error code: ', blacs_info
      call messages_fatal(1, namespace=st%namespace)
    end if
     
    wsize = nint(R_REAL(tmp))
    SAFE_ALLOCATE(work(1:wsize))
#ifdef HAVE_SCALAPACK
#ifndef R_TREAL
    rwsize = max(1,nint(R_REAL(tmp2)))
    SAFE_ALLOCATE(rwork(1:rwsize))
    call pzgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), work(1), wsize, rwork(1), rwsize, blacs_info)
    SAFE_DEALLOCATE_A(rwork)
#else
    call pdgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), work(1), wsize,  blacs_info)
#endif
#endif

    if(blacs_info /= 0) then
      write(message(1),'(a,i6)') 'scalapack geqrf call failed with error code: ', blacs_info
      call messages_fatal(1, namespace=st%namespace)
    end if
    SAFE_DEALLOCATE_A(work)
     
     ! copy the first nst global elements of ipiv into jpvt
     ! bcast is at the end of the routine
!     if(mpi_world%rank==0)  then
!        do ist =1,nst
!           write(123,*) ipiv(ist)
!        end do
!     end if
    jpvt(1:nst) =  ipiv(1:nst)
     
  else
    ! first gather states into one array on the root process
    ! build transpose of KS set on which RRQR is performed
    if(root) then
       SAFE_ALLOCATE(KSt(1:nst,1:mesh%np_global))
    end if
    
    ! gather states in case of domain parallelization
    if (mesh%parallel_in_domains.or.st%parallel_in_states) then
      SAFE_ALLOCATE(temp_state(1:mesh%np,1))
      SAFE_ALLOCATE(state_global(1:mesh%np_global))
      
      count = 0
      do ii = 1,nst
        !we are copying states like this:  KSt(i,:) = st%psi(:,dim,i,nik)
        state_global(1:mesh%np_global) = M_ZERO
#ifdef HAVE_MPI
        sender = 0
        if(state_is_local(st,ii)) then
          call states_elec_get_state(st, mesh, ii, ik, temp_state)
          call vec_gather(mesh%vp, 0, temp_state(1:mesh%np,1), state_global)
          if(mesh%mpi_grp%rank ==0) sender = mpi_world%rank
        end if
        call comm_allreduce(mpi_world%comm,sender)
        call MPI_Bcast(state_global,mesh%np_global , R_MPITYPE, sender, mpi_world%comm, mpi_err)
#endif
        ! keep full Kohn-Sham matrix only on root
        if (root)  KSt(ii,1:mesh%np_global)  = st%occ(ii,1)*state_global(1:mesh%np_global)
      end do
      SAFE_DEALLOCATE_A(state_global)
      SAFE_DEALLOCATE_A(temp_state)
    else
      ! serial
      SAFE_ALLOCATE(temp_state(1:mesh%np,1))
      do ii = 1, nst
        ! this call is necessary becasue we want to have only np not np_part
        call states_elec_get_state(st, mesh, ii, ik, temp_state)
        KSt(ii,:) = st%occ(ii,1)*temp_state(:,1)
      end do
      SAFE_DEALLOCATE_A(temp_state)
    end if

    ! now perform serial RRQR
    ! dummy call to obtain dimension of work
    ! Note: the lapack routine has different number of arguments depending on type
    if(root) then
      SAFE_ALLOCATE(work(1:1))
      SAFE_ALLOCATE(tau(1:nst))
      if(.not.states_are_real(st)) then
         SAFE_ALLOCATE(rwork(1:2*mesh%np_global))
         call zgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, -1, rwork, info)
      else
         call dgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, -1, info)
      endif
      if (info /= 0) then
         write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
         call messages_fatal(1, namespace=st%namespace)
      end if

      wsize = work(1)
      SAFE_DEALLOCATE_A(work)
      SAFE_ALLOCATE(work(1:wsize))

      jpvt(:) = 0
      tau(:) = 0.
      ! actual call
      if(.not.states_are_real(st)) then
         call zgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, wsize, rwork, info)
      else
         call dgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, wsize, info)
      endif
      if (info /= 0)then
         write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
         call messages_fatal(1, namespace=st%namespace)
      end if
      SAFE_DEALLOCATE_A(work)
    endif

    SAFE_DEALLOCATE_A(temp_state)
    SAFE_DEALLOCATE_A(state_global)
    
   endif

#ifdef HAVE_MPI
    call MPI_Bcast(JPVT,nst, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
#endif

   call profiling_out(prof_scdm_QR)
   POP_SUB(X(scdm_rrqr))

end subroutine X(scdm_rrqr)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
