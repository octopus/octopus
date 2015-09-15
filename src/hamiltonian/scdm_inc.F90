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
!! $Id$

!> this performs the SCDM localization, transforming the original set of states KS (Kohn-Sham)
!! into the set SCDM, by first performing RRQR and then Cholesky for orthogonalization
subroutine X(scdm_localize)(st,mesh,scdm)

  type(states_t), intent(in) :: st !< this contains the non-localize set KS (for now from hm%hf_st which is confusing)
  type(mesh_t), intent(in)   :: mesh
  type(scdm_t) :: scdm

  integer :: ii, jj, kk, ll, vv, count, ip, nval, info,i1, i2, i3, idim, j1, j2, j3
  integer :: JPVT(mesh%np_global)
  integer :: icenter(3), ind_center

  integer :: nn(3)
  R_TYPE, allocatable :: KSt(:,:), KSt_original(:,:)
  R_TYPE, allocatable :: SCDM_temp(:,:), Pcc(:,:)
  R_TYPE, allocatable :: rho(:), pot(:), rho2(:)
  FLOAT  :: exx, error, error_tmp
  R_TYPE, allocatable :: state_global(:), temp_state(:,:)
  integer,  allocatable :: temp_box(:,:,:)
  integer :: ix(3), nr(3,2)
  logical :: out_of_index_range(3), out_of_mesh(3)
  FLOAT, allocatable :: lxyz_domains(:,:),  lxyz_global(:), lxyz_local(:)
  
  FLOAT :: t0,t1, t2, t3, t4
  FLOAT :: temp(3)
  character(len=50) :: name
  type(cube_function_t) :: cf

  PUSH_SUB(X(scdm_localize))
  ! check if already localized
  if(scdm_is_local) then
    POP_SUB(X(scdm_localize))
    return
  end if

  call cpu_time(t0)
  if (st%lnst /= st%nst) call messages_not_implemented("SCDM with state parallelization")
  nval = st%nst ! TODO: check that this is really the number of valence states

  ! build transpose of KS set on which RRQR is performed
  if (scdm%root) then
    SAFE_ALLOCATE(KSt(1:nval,1:mesh%np_global))
    ! keep a copy of this NOTE: maybe too expensive in memory?
    SAFE_ALLOCATE(KSt_original(1:nval,1:mesh%np_global))
  else
    ! non root processes hold their part of the Kohn-Sham states for distributed sum
    SAFE_ALLOCATE(KSt(scdm%lnst,mesh%np_global))
    SAFE_ALLOCATE(KSt_original(scdm%lnst,mesh%np_global))
  end if

  !NOTE: not sure how to proceed if dim!=1 or nik!=1
  if (st%d%nik /= 1 .or. st%d%dim /= 1) call messages_not_implemented("SCDM with k-points or dims")
  ! gather states in case of domain parallelization
  if (mesh%parallel_in_domains) then
    SAFE_ALLOCATE(temp_state(1:mesh%np,1))
    SAFE_ALLOCATE(state_global(1:mesh%np_global))

    count = 0
    do ii = 1, nval
      ! KSt(i,:) = st%psi(:,dim,i,nik)
#ifdef HAVE_MPI
      call states_get_state(st, mesh, ii, st%d%nik, temp_state)
      call vec_gather(mesh%vp, 0, state_global, temp_state(1:mesh%np,1))  
      call MPI_Bcast(state_global,mesh%np_global , R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif
      if (scdm%root) then
        KSt(ii,:)  = st%occ(ii,1)*state_global(:)
      else
        ! on non-root processes keep only local number of Kohn-Sham states
        if (ii >= scdm%st_start .and. ii <= scdm%st_end) then
          count = count + 1
          KSt(count,1:mesh%np_global) = st%occ(ii,1)*state_global(1:mesh%np_global)
        end if
      end if
    end do
    SAFE_DEALLOCATE_A(state_global)
    SAFE_DEALLOCATE_A(temp_state)
  else
    ! serial
    SAFE_ALLOCATE(temp_state(1:mesh%np,1))
    do ii = 1, nval
      ! this call is necessary becasue we want to have only np not np_part
      call states_get_state(st, mesh, ii, st%d%nik, temp_state)
      KSt(ii,:) = st%occ(ii,1)*temp_state(:,1)
    end do
    SAFE_DEALLOCATE_A(temp_state)
  end if

  ! possibly redundant copy
  KSt_original(:,:) = KSt(:,:)

  call cpu_time(t1)

  if(scdm%root) then

    call X(RRQR)(nval,mesh%np_global,KSt,JPVT)
    call cpu_time(t2)
    if(scdm%verbose) call messages_print_var_value(stdout, 'time: RRQR:', t2-t1)
    
  end if
#ifdef HAVE_MPI
  call MPI_Bcast(JPVT,mesh%np_global , MPI_INTEGER, 0, mesh%mpi_grp%comm, mpi_err)
#endif

  SAFE_DEALLOCATE_A(KSt)

  ! form SCDM, Note: This could be done in one step together with the orthogonalization
  !                  to save this allocation
  SAFE_ALLOCATE(SCDM_temp(1:mesh%np_global,1:nval))
  SCDM_temp(:,:) = M_ZERO
  do ii = 1, nval
    do vv = 1, scdm%lnst
      SCDM_temp(1:mesh%np_global,ii) = SCDM_temp(1:mesh%np_global,ii) + &
                 KSt_original(vv,1:mesh%np_global)*R_CONJ( KSt_original(vv,JPVT(ii)) )
    end do
  end do
  
  call comm_allreduce(mesh%mpi_grp%comm, SCDM_temp, (/mesh%np_global, nval/))

  call cpu_time(t1)
  if (scdm%verbose) call messages_print_var_value(stdout, 'time: explicit matmul1:',t1-t2)

    ! --- Orthogoalization ----
    ! form lower triangle of Pcc
    SAFE_ALLOCATE(Pcc(1:nval,1:nval))
    Pcc(:,:) = M_ZERO
    ! work only on root process, because its the only one that holds the full KSt_original
    if(scdm%root) then
      do ii = 1, nval
        do jj = 1, ii
          do vv = 1, nval
            Pcc(ii,jj) = Pcc(ii,jj)+ KSt_original(vv,JPVT(ii))*R_CONJ(KSt_original(vv,JPVT(jj)))
          end do
        end do
      end do

      call cpu_time(t2)
      if(scdm%verbose) call messages_print_var_value(stdout, 'time: explicit matmul2:',t2-t1)
      ! Cholesky fact.
      call X(POTRF)("L", nval, Pcc, nval, info )
      if (info /= 0) then
        if (info < 0) then
          write(message(1),'(A28,I2)') 'Illegal argument in DPOTRF: ', info
          call messages_fatal(1)
        else
          message(1) = 'Fail of Cholesky, not pos-semi-def '
          call messages_fatal(1)
        end if
        stop
      end if

      call cpu_time(t1)
      if(scdm%verbose) call messages_print_var_value(stdout, 'time: cholesky:',t1-t2)
      ! transpose
      Pcc(:,:) = transpose(R_CONJ(Pcc(:,:)))
      ! invert
      call X(invert)(nval,Pcc)
    end if
#ifdef HAVE_MPI
    call MPI_Bcast(Pcc,nval*nval , R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif
    call cpu_time(t2)
    if(scdm%verbose) call messages_print_var_value(stdout, 'time: transpose invert:',t2-t1)

    ! form orthogonal SCDM
    ! NOTE: this needs state parallelization
    SAFE_ALLOCATE(temp_state(1:mesh%np,1))
    SAFE_ALLOCATE(state_global(1:mesh%np_global))
    do vv = 1, nval ! st%start ... st%end do only local state distributed part of the product
      state_global(:) = M_ZERO
      do jj = 1, nval
        state_global(1:mesh%np_global) = state_global(1:mesh%np_global) + SCDM_temp(1:mesh%np_global,jj)*Pcc(jj, vv)
      end do
      call vec_scatter(mesh%vp, 0, state_global, temp_state(1:mesh%np,1))
      call states_set_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
    end do
    
    SAFE_DEALLOCATE_A(SCDM_temp)
    
    call cpu_time(t1)
    if (scdm%verbose) call messages_print_var_value(stdout,  'time: explicit matmul3',t1-t2)

    ! normalise SCDM states
    do vv = 1, nval!scdm%lnst ! this needs state parallelization
      call states_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
      temp_state(1:mesh%np,:) = temp_state(1:mesh%np,:)/X(mf_nrm2)(mesh,temp_state(1:mesh%np,1))
      call states_set_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
    end do
    call cpu_time(t2)
    if(scdm%verbose) call messages_print_var_value(stdout,  'time: norms',t2-t1)

    ! check orthonormality
    !print *, 'orthonrmality: ================'
    !do j=1,nval
    !   do i=1,nval
    !      print *, i,j, &
    !        dot_product(scdm%st%X(dontusepsi)(1:mesh%np,1,i,1),scdm%st%X(dontusepsi)(1:mesh%np,1,j,1))*mesh%volume_element
    !      print *, i,j, dot_product(st%X(dontusepsi)(1:mesh%np,1,i,1),st%X(dontusepsi)(1:mesh%np,1,j,1))*mesh%volume_element
    !   end do
    !end do
    !print *, '=============================='

    ! write cube files
    !call X(io_function_output)(io_function_fill_how('Cube'), ".", "SCDM_1", mesh, scdm%st%X(dontusepsi)(:,1,1,1), &
    !                            unit_one, info,geo=scdm_geo)

    call cpu_time(t1)
    !       print *, 'time: output',t1-t2

    ! end of SCDM procedure
    
    ! find centers, by computing center of mass of |psi|^2
    scdm%center(:,:) = 0
    count = 0
    ! get domains of mesh%idx%lxyz
    SAFE_ALLOCATE(lxyz_domains(1:mesh%np,3))
    SAFE_ALLOCATE(lxyz_global(1:mesh%np_global))
    SAFE_ALLOCATE(lxyz_local(1:mesh%np))
    do ii=1,3
      lxyz_global(1:mesh%np_global) = mesh%idx%lxyz(1:mesh%np_global,ii)
      call vec_scatter(mesh%vp, 0, lxyz_global, lxyz_local)
      lxyz_domains(1:mesh%np,ii) = lxyz_local(1:mesh%np)
    end do
    
    do vv = 1, nval
      call states_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
      do ii=1,3
        scdm%center(ii,vv) = sum(temp_state(1:mesh%np,1)*R_CONJ(temp_state(1:mesh%np,1))*&
                          lxyz_domains(1:mesh%np,ii)*mesh%spacing(ii))*mesh%volume_element
      end do
    end do
    call comm_allreduce(mesh%mpi_grp%comm, scdm%center)

    SAFE_DEALLOCATE_A(lxyz_domains)
    SAFE_DEALLOCATE_A(lxyz_global)
    SAFE_DEALLOCATE_A(lxyz_local)
    call cpu_time(t2)
    if (scdm%verbose) call messages_print_var_value(stdout, 'time: find centers',t2-t1)

    ! end of SCDM procedure
#ifdef HAVE_MPI
  call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif

  ! copy local box of state
  call cpu_time(t1)
  count = 0
  error_tmp = M_ZERO
  scdm%X(psi)(:,:) =  M_ZERO
  SAFE_ALLOCATE(temp_box(1:2*scdm%box_size+1,1:2*scdm%box_size+1,1:2*scdm%box_size+1))
  do vv = 1,nval! Needs state distribution: scdm%st_start, scdm%st_end
    count = count + 1
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
        call messages_fatal(1)
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
            call messages_fatal(1)
          end if
        end if
      end do
      
    end if
    !
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
      !print *, vv, minval(temp_box), maxval(temp_box),  mesh%np_global
      !do i1=1,scdm%box_size*2+1
      !  do i2=1,scdm%box_size*2+1
      !    do i3=1,scdm%box_size*2+1
      !      write(333,*) real(temp_box(i1,i2,i3))
      !    end do
      !  end do
      !end do
      message(1) = 'SCDM box mapping failed'
      call messages_fatal(1)
    end if
    
    scdm%box(:,:,:,count) = temp_box(:,:,:)

    ! to copy the scdm state in the box, we need the global state
    call states_get_state(scdm%st, mesh, vv, scdm%st%d%nik, temp_state(1:mesh%np,:))
    call vec_gather(mesh%vp, 0, state_global, temp_state(1:mesh%np,1))
#ifdef HAVE_MPI
    call MPI_Bcast(state_global,mesh%np_global , R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif    
    ! copy points to box
    ! this box refers to the global mesh
    do jj = 1, scdm%box_size*2+1
      do kk = 1, scdm%box_size*2+1
        do ll = 1, scdm%box_size*2+1
          ! map into the box
          ip = (jj-1)*((scdm%box_size*2+1))**2+(kk-1)*((scdm%box_size*2+1)) + ll
          scdm%X(psi)(ip,count) = state_global(scdm%box(jj,kk,ll,count))
        end do
      end do
    end do

    ! compute localization error
    error_tmp = error_tmp + M_ONE - dot_product(scdm%X(psi)(:,count),scdm%X(psi)(:,count))*mesh%volume_element

    ! re-normalize inside box
    ! if(scdm%re_ortho_normalize) then
    !    scdm%X(psi)(:,count) = scdm%X(psi)(:,count)/(dot_product(scdm%X(psi)(:,count),scdm%X(psi)(:,count))*mesh%volume_element)
    !    !
    !    ! for testing zero outside the box
    !    scdm%st%X(dontusepsi)(:,st%d%dim,vv,scdm%st%d%nik) = 0.
    !    do jj=1,scdm%box_size*2+1
    !       do kk=1,scdm%box_size*2+1
    !          do ll=1,scdm%box_size*2+1
    !             ip = (jj-1)*((scdm%box_size*2+1))**2+(kk-1)*((scdm%box_size*2+1)) + ll
    !             scdm%st%X(dontusepsi)(scdm%box(jj,kk,ll,vv),st%d%dim,vv,scdm%st%d%nik) = scdm%X(psi)(ip,count)
    !          end do
    !       end do
    !    end do
    ! end if
    
  end do

  SAFE_DEALLOCATE_A(temp_box)
  
#ifdef HAVE_MPI
  error = M_ZERO
  call MPI_Allreduce(error_tmp, error, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
#else
  error = error_tmp
#endif
  if (scdm%root .and. scdm%verbose) call messages_print_var_value(stdout, 'SCDM localization error:', error/st%nst)
  !
  call cpu_time(t2)
  if (scdm%root .and. scdm%verbose) call messages_print_var_value(stdout, 'time: copy box',t2-t1)

  !if(scdm%re_ortho_normalize) then
  !   ! check orthonormality
  !   print *, 'orthonrmality in boxes: ================'
  !   do j=1,nval
  !      do i=1,nval
  !         print *, i,j, &
  !           dot_product(scdm%st%X(dontusepsi)(1:mesh%np,1,i,1),scdm%st%X(dontusepsi)(1:mesh%np,1,j,1))*mesh%volume_element
  !      end do
  !   end do
  !   print *, '=============================='
  !   ! and re-orthogonalize-----------------------------------------------------
  !   call X(states_orthogonalization_full)(scdm%st, mesh, 1)
  !   ! and cut again
  !   count = 0
  !   scdm%X(psi)(:,:) =  M_ZERO
  !   do v=scdm%st_start,scdm%st_end
  !      count = count +1
  !       do i=1,3
  !          icenter(i) = scdm%center(i,v)/mesh%spacing(i)
  !       end do
  !       ind_center = mesh%idx%lxyz_inv(icenter(1),icenter(2),icenter(3))
  !       scdm%box(:,:,:,count) =  mesh%idx%lxyz_inv(icenter(1)-scdm%box_size:icenter(1)+scdm%box_size, &
  !                                              icenter(2)-scdm%box_size:icenter(2)+scdm%box_size, &
  !                                              icenter(3)-scdm%box_size:icenter(3)+scdm%box_size)
  !       do j=1,scdm%box_size*2+1
  !          do k=1,scdm%box_size*2+1
  !             do l=1,scdm%box_size*2+1
  !                ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
  !                scdm%X(psi)(ip,count) = scdm%st%X(dontusepsi)(scdm%box(j,k,l,count),st%d%dim,count,scdm%st%d%nik)
  !             end do
  !          end do
  !       end do
  !       ! and set to zero outside again
  !       scdm%st%X(dontusepsi)(:,st%d%dim,v,scdm%st%d%nik) = 0.
  !       do j=1,scdm%box_size*2+1
  !          do k=1,scdm%box_size*2+1
  !             do l=1,scdm%box_size*2+1
  !                ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
  !                scdm%st%X(dontusepsi)(scdm%box(j,k,l,v),st%d%dim,v,scdm%st%d%nik) = scdm%X(psi)(ip,count)
  !             end do
  !          end do
  !       end do
  !       !
  !    end do
  !!--------------------------------------------------------------------------
  !print *, 'orthonrmality in boxes after re-ortho: ================'
  !do j=1,nval
  !   do i=1,nval
  !      print*, i,j, &
  !        dot_product(scdm%st%X(dontusepsi)(1:mesh%np,1,i,1),scdm%st%X(dontusepsi)(1:mesh%np,1,j,1))*mesh%volume_element
  !   end do
  !end do
  !print *, '======================================================'
  !end if
  !

  ! check span
  !SAFEx_ALLOCATE(state_global(1:mesh%np_global))
  !state_global(:) = M_ZERO
  !print *, 'quality of projector:'
  !do j=1,nval
  !   state_global(:) = M_ZERO
  !   do i=1,nval
  !      state_global(:) = state_global(:) + &
  !           dot_product(scdm%st%X(dontusepsi)(:,1,i,1),st%X(dontusepsi)(1:mesh%np,1,j,1))&
  !           *scdm%st%X(dontusepsi)(:,1,i,1)*mesh%volume_element
  !   end do
  !   print *, j, M_ONE - X(mf_nrm2)(mesh, state_global)
  !end do

  ! set flag to do this only once
  scdm_is_local = .true.

  SAFE_DEALLOCATE_A(Pcc)
  SAFE_DEALLOCATE_A(KSt)
  SAFE_DEALLOCATE_A(scdm_temp)

  call cpu_time(t1)
  if (scdm%root .and. scdm%verbose) call messages_print_var_value(stdout,  'time: all SCDM',t1-t0)

  !    return

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  ! this is to comupte exchange energy, in priciple... not operational
  !   SxAFE_ALLOCATE(rho(1:mesh%np_part))
  !   SxAFE_ALLOCATE(rho2(1:mesh%np_part))
  !   SxAFE_ALLOCATE(pot(1:mesh%np_part))
  !   !
  !   exx = 0.
  !   do i=1,nval
  !      do j=1,i
  !         temp(:) = (scdm%center(:,i)- scdm%center(:,j))/mesh%spacing(:)
  !!          if(sqrt(dot_product(temp,temp)).le.2.*scdm%rcut) then
  !            !
  !            rho(:) = scdm%st%X(dontusepsi)(:,st%d%dim,i,scdm%st%d%nik)*scdm%st%X(dontusepsi)(:,st%d%dim,j,scdm%st%d%nik)
  !!rho(:) = st%X(dontusepsi)(:,st%d%dim,i,scdm%st%d%nik)*st%X(dontusepsi)(:,st%d%dim,j,scdm%st%d%nik)
  !            rho2(:) = rho(:)
  !            pot(:) = 0.
  !            call X(poisson_solve)(psolver, pot, rho, all_nodes = .false.)
  !            !
  !            ! catch diagonals for double counting
  !            if(i.ne.j) then
  !               exx = exx - 0.5*dot_product(pot(:),rho2(:))*mesh%volume_element
  !            else
  !               exx = exx - 0.25*dot_product(pot(:),rho2(:))*mesh%volume_element
  !            end if
  !!          end if
  !         !
  !      end do
  !   end do
  !   !
  !   call messages_print_var_value(stdout,'exx[eV] = ', exx*27.211396132)
  !   SxAFE_DEALLOCATE_A(rho)
  !   SxAFE_DEALLOCATE_A(rho2)
  !   SxAFE_DEALLOCATE_A(pot)

  POP_SUB(X(scdm_localize))

end subroutine X(scdm_localize)

!> rotate states from KS to SCDM representation and construct SCDM states
subroutine X(scdm_rotate_states)(st,mesh,scdm)
  type(states_t), intent(inout)  :: st
  type(mesh_t), intent(in)       :: mesh
  type(scdm_t), intent(inout)    :: scdm

  PUSH_SUB(X(scdm_rotate_states))

  ! create localized SCDM representation of the states in st
  scdm_is_local = .false.
  call X(scdm_localize)(st,mesh,scdm)

  ! overwrite state object with the scdm states
  call states_copy(st,scdm%st)

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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
