  subroutine X(scdm_localize)(st,mesh,scdm)
    !
    ! this performs the SCDM localization, transforming the original set of states KS (Kohn-Sham)
    ! into the set SCDM, by first performing RRQR and then Cholesky for orthogonalization
    !
    type(states_t), intent(in)   :: st ! this contains the non-localize set KS (for now from hm%hf_st which is confusing)
    type(mesh_t), intent(in)     :: mesh
    type(scdm_t) :: scdm
    !
    integer            :: i,j,k,l,v,count,ip, nval, INFO,i1, i2, i3, idim, j1,j2,j3
    integer            ::JPVT(mesh%np_global)
    integer            :: icenter(3), ind_center
    !
integer :: nn(3)
    R_TYPE, allocatable :: KSt(:,:), KSt_original(:,:)
    R_TYPE, allocatable ::  SCDM_temp(:,:), Pcc(:,:)
    R_TYPE, allocatable :: rho(:), pot(:), rho2(:)
    FLOAT  :: exx, error, error_tmp
    R_TYPE, allocatable     :: state_global(:), temp_state(:,:)
    integer :: ix(3)
    logical  :: outside
    !
    FLOAT :: t0,t1, t2, t3, t4
    FLOAT :: temp(3)
    character(len=50) :: name
    type(cube_function_t) :: cf
    !
    PUSH_SUB(X(scdm_localize))
    ! check if already localized
    if(scdm_is_local) return
    !
    call cpu_time(t0)
    if(st%lnst.ne.st%nst) call messages_not_implemented("SCDM with state parallelization")
    nval = st%nst ! TODO: check that this is really the number of valence states
    !
    ! built transpose of KS set on which RRQR is performed
    if(scdm%root) then
       SAFE_ALLOCATE(KSt(nval,mesh%np_global))
       ! keep a copy of this NOTE: maybe too expensive in memory?
       SAFE_ALLOCATE(KSt_original(nval,mesh%np_global))
    endif
    !
    !NOTE: not sure how to proceed if dim!=1 or nik!=1
    if( st%d%nik.ne.1.or.st%d%dim.ne.1) call messages_not_implemented("SCDM with k-points or dims")
    ! gather states in case of domain paralleization
    if(mesh%parallel_in_domains) then
       SAFE_ALLOCATE(state_global(1:mesh%np_global))
       !
       do i=1,nval
          ! KSt(i,:) = st%dpsi(:,st%d%dim,i,st%d%nik)
#if HAVE_MPI
          call vec_gather(mesh%vp, 0, state_global, st%X(psi)(1:mesh%np,st%d%dim,i,st%d%nik))
#endif
          if(scdm%root) KSt(i,:)  = st%occ(i,1)*state_global(:)
       end do
       SAFE_DEALLOCATE_A(state_global)
    else
       ! serial
       SAFE_ALLOCATE(temp_state(1:mesh%np,1))
       do i=1,nval
          ! this call is necessary becasue we want to have only np not np_part
           call states_get_state(st, mesh, i, st%d%nik, temp_state)
           KSt(i,:) = st%occ(i,1)*temp_state(:,1)
       end do
       SAFE_DEALLOCATE_A(temp_state)
    endif
    !
    ! possibly redundant copy
    if(scdm%root) KSt_original(:,:) = KSt(:,:)
    !
    !
    call cpu_time(t1)
    ! perform the RRQR
    scdm%st%X(psi)(:,:,:,:) = M_ZERO ! this is important for distribution later
    ! Note: this should be parallelized, once it is clear that we would save time 
!---! ----------------------------- SERIAL START --------------------------------
    if(scdm%root) then
       !
       call X(RRQR)(nval,mesh%np_global,KSt,JPVT)
       call cpu_time(t2)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: RRQR:', t2-t1)
       !
       SAFE_DEALLOCATE_A(KSt)
       !
       ! form SCDM Note: This could be done in one step together with the orhtogonalization
       !                 to save this allocation
       SAFE_ALLOCATE(SCDM_temp(mesh%np_global,nval))
       SCDM_temp(:,:) = M_ZERO
       do i=1,nval
          do v=1,nval
             SCDM_temp(:,i) = SCDM_temp(:,i) + KSt_original(v,:)*R_CONJ(KSt_original(v,JPVT(i)))
          end do
       end do
       !
       call cpu_time(t1)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: explicit matmul1:',t1-t2)
       !
       ! --- Orthogoalization ----
       ! form lower trinagle of Pcc
       SAFE_ALLOCATE(Pcc(nval,nval))
       Pcc(:,:) = M_ZERO
       do i=1,nval
          do j=1,i
             do v=1,nval
                Pcc(i,j) = Pcc(i,j)+ KSt_original(v,JPVT(i))*R_CONJ(KSt_original(v,JPVT(j)))
             enddo
          enddo
       end do
       !
       call cpu_time(t2)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: explicit matmul2:',t2-t1)
       ! Cholesky fact.
       call X(POTRF)("L", nval, Pcc, nval, INFO )
       if(INFO.ne.0) then
          if(INFO.lt.0) then
             print *, 'Illegal argument in DPOTRF: ', INFO
          else
             print *, 'Fail of Cholesky, not pos-semi-def '
          endif
          stop
       endif
       !
       call cpu_time(t1)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: cholesky:',t1-t2)
       ! transpose
       Pcc(:,:) = transpose(R_CONJ(Pcc(:,:)))
       ! invert
       call X(invert)(nval,Pcc)
       !
       call cpu_time(t2)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: transpose invert:',t2-t1)
       ! form orthogonal SCDM
       scdm%st%X(psi)(:,:,:,:) = M_ZERO
       do i=1,mesh%np_global
          do j=1,nval
             do v=1,nval
                scdm%st%X(psi)(i,1,v,1) = scdm%st%X(psi)(i,1,v,1) + SCDM_temp(i,j)*Pcc(j,v)
             end do
          end do
       end do
       !
       SAFE_DEALLOCATE_A(SCDM_temp)
       call cpu_time(t1)
       if(scdm%verbose) call messages_print_var_value(stdout,  'time: explicit matmul3',t1-t2)
       ! normalise SCDM states
       do v=1,nval
          scdm%st%X(psi)(:,1,v,1) = scdm%st%X(psi)(:,1,v,1)/&
               (sqrt(dot_product(scdm%st%X(psi)(:,1,v,1),scdm%st%X(psi)(:,1,v,1))*mesh%volume_element))
               !this should be used ../X(mf_nrm2)(mesh,scdm%st%X(psi)(:,1,v,1))
               ! but doensnt work with parallelization
       end do
       call cpu_time(t2)
       if(scdm%verbose) call messages_print_var_value(stdout,  'time: norms',t2-t1)
       !
       ! check orthonormality
       !print *, 'orthonrmality: ================'
       !do j=1,nval
       !   do i=1,nval
       !      print *, i,j, dot_product(scdm%st%X(psi)(1:mesh%np,1,i,1),scdm%st%X(psi)(1:mesh%np,1,j,1))*mesh%volume_element
       !      print *, i,j, dot_product(st%X(psi)(1:mesh%np,1,i,1),st%X(psi)(1:mesh%np,1,j,1))*mesh%volume_element
       !   enddo
       !enddo
       !print *, '=============================='
       !
       ! write cube files
       !call X(io_function_output)(io_function_fill_how('Cube'), ".", "SCDM_1", mesh, scdm%st%X(psi)(:,1,1,1), &
       !                            unit_one, info,geo=scdm_geo)
       !
       call cpu_time(t1)
!       print *, 'time: output',t1-t2
       !
       ! find centers, by computing center of mass of |psi|^2
       scdm%center(:,:) = 0
       do v=1,nval
          do i=1,3
             scdm%center(i,v) = sum(scdm%st%X(psi)(:,st%d%dim,v,scdm%st%d%nik)*R_CONJ(scdm%st%X(psi)(:,st%d%dim,v,scdm%st%d%nik))* &
                       mesh%idx%lxyz(1:mesh%np_global,i)*mesh%spacing(i))*mesh%volume_element
          end do
          write(127,*) scdm%center(:,v)
       end do
       close(127)
       call cpu_time(t2)
       if(scdm%verbose) call messages_print_var_value(stdout, 'time: find centers',t2-t1)
    endif
    !
!---! --------------------- SERIAL END ------------------------------------
    !
#if HAVE_MPI
    call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif
    !
    ! distribute the localized states
    count = 0
    SAFE_ALLOCATE(temp_state(1:mesh%np_global,1))
    do i=1,nval
       ! send state to all processes as temp
       temp_state(:,:) = M_ZERO
       ! this is only non-zero on root
       if(scdm%root) temp_state(:,1) =  scdm%st%X(psi)(1:mesh%np_global,st%d%dim,i,scdm%st%d%nik)
#if HAVE_MPI
       call MPI_Bcast(temp_state(1,1),mesh%np_global , R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif
       !
       ! only keep the state if it falls into index range on process
       if(i.ge.scdm%st_start.and.i.le.scdm%st_end) then
          count = count +1
          scdm%st%X(psi)(1:mesh%np_global,st%d%dim,count,scdm%st%d%nik) = temp_state(:,1)
       end if
       !
    end do
    SAFE_DEALLOCATE_A(temp_state)
    !
#if HAVE_MPI
    ! broadcast the centers to all processes
    call MPI_Bcast(scdm%center(1,1),size(scdm%center) ,MPI_FLOAT, 0, mesh%mpi_grp%comm, mpi_err)
#endif
    !
    ! copy local box of state
    call cpu_time(t1)
    count = 0
    error_tmp = M_ZERO
    scdm%X(psi)(:,:) =  M_ZERO
    do v=scdm%st_start,scdm%st_end
       count = count +1
       ! find integer index of center
       do i=1,3
          icenter(i) = scdm%center(i,v)/mesh%spacing(i)
       end do
       ! find index of center in the mesh
       ind_center = mesh%idx%lxyz_inv(icenter(1),icenter(2),icenter(3))
       !
       ! if(periodic boundary..)
       call check_periodic_box(mesh%idx,icenter(:),scdm%box_size,scdm%periodic(count))
       !
       ! make list with points in the box
       if(.not.scdm%periodic(count)) then
          scdm%box(:,:,:,count) =  mesh%idx%lxyz_inv(icenter(1)-scdm%box_size:icenter(1)+scdm%box_size, &
                                              icenter(2)-scdm%box_size:icenter(2)+scdm%box_size, &
                                              icenter(3)-scdm%box_size:icenter(3)+scdm%box_size)
       else
          print *, 'copy periodic state', count
          ! in case there are periodic replica go through every point
          do i1=-scdm%box_size,scdm%box_size
             do i2=-scdm%box_size,scdm%box_size
                do i3=-scdm%box_size,scdm%box_size
                   !
                   ix(:) = icenter(:)+(/i1,i2,i3/)
                   !
                   outside = .false.
                   do idim=1,3
                      if(ix(idim).lt.mesh%idx%nr(1,idim).or.ix(idim).gt.mesh%idx%nr(2,idim)) then
                         outside = .true.
                         exit
                      end if
                   end do
                   !
                   if(outside) then 
                      ! map point to equivalent one in cell
                      do idim=1,3 ! should be dim...
                         if( ix(idim).lt.mesh%idx%nr(1,idim) ) then
                            ix(idim) = ix(idim)+mesh%idx%ll(idim)
                         elseif( ix(idim).gt.mesh%idx%nr(2,idim) ) then
                            ix(idim) = ix(idim)-mesh%idx%ll(idim)
                         end if
                      end do
                   end if
                   ! indices of box are 1-based
                   j1=i1+scdm%box_size+1
                   j2=i2+scdm%box_size+1
                   j3=i3+scdm%box_size+1
                   scdm%box(j1,j2,j3,count) = mesh%idx%lxyz_inv(ix(1),ix(2),ix(3))
                   !
                end do!i3
             end do!i2
          end do!i1
          !
       end if
       !
       ! copy points to box
       ! this box refers to the global mesh
       do j=1,scdm%box_size*2+1
          do k=1,scdm%box_size*2+1
             do l=1,scdm%box_size*2+1
                ! map into the twice larger box
!                ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
ip = (j-1)*((scdm%box_size*2+1))**2+(k-1)*((scdm%box_size*2+1)) + l
                scdm%X(psi)(ip,count) = scdm%st%X(psi)(scdm%box(j,k,l,count),st%d%dim,count,scdm%st%d%nik)
             end do
          end do
       end do
       !
       ! compue localization error
       error_tmp = error_tmp + M_ONE - dot_product(scdm%X(psi)(:,count),scdm%X(psi)(:,count))*mesh%volume_element
       !
       ! re-normalize inside box
!       if(scdm%re_ortho_normalize) then
!          scdm%X(psi)(:,count) = scdm%X(psi)(:,count)/(dot_product(scdm%X(psi)(:,count),scdm%X(psi)(:,count))*mesh%volume_element)
!          !
!          ! for testing zero outside the box
!          scdm%st%X(psi)(:,st%d%dim,v,scdm%st%d%nik) = 0.
!          do j=1,scdm%box_size*2+1
!             do k=1,scdm%box_size*2+1
!                do l=1,scdm%box_size*2+1
!                   ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
!                   scdm%st%X(psi)(scdm%box(j,k,l,v),st%d%dim,v,scdm%st%d%nik) = scdm%X(psi)(ip,count)
!                enddo
!             enddo
!          enddo
!       endif
       !
   end do
   !
#if HAVE_MPI
   error = M_ZERO
   call MPI_Allreduce(error_tmp, error, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
#else
   error = error_tmp
#endif
   if(scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'SCDM localization error:', error/st%nst)
   !
   call cpu_time(t2)
   if(scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'time: copy box',t2-t1)
   !
!
!if(scdm%re_ortho_normalize) then
!   ! check orthonormality
!   print *, 'orthonrmality in boxes: ================'
!   do j=1,nval
!      do i=1,nval
!         print *, i,j, dot_product(scdm%st%X(psi)(1:mesh%np,1,i,1),scdm%st%X(psi)(1:mesh%np,1,j,1))*mesh%volume_element
!      enddo
!   enddo
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
!       enddo
!       ind_center = mesh%idx%lxyz_inv(icenter(1),icenter(2),icenter(3))
!       scdm%box(:,:,:,count) =  mesh%idx%lxyz_inv(icenter(1)-scdm%box_size:icenter(1)+scdm%box_size, &
!                                              icenter(2)-scdm%box_size:icenter(2)+scdm%box_size, &
!                                              icenter(3)-scdm%box_size:icenter(3)+scdm%box_size)
!       do j=1,scdm%box_size*2+1
!          do k=1,scdm%box_size*2+1
!             do l=1,scdm%box_size*2+1
!                ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
!                scdm%X(psi)(ip,count) = scdm%st%X(psi)(scdm%box(j,k,l,count),st%d%dim,count,scdm%st%d%nik)
!             enddo
!          enddo
!       enddo
!       ! and set to zero outside again
!       scdm%st%X(psi)(:,st%d%dim,v,scdm%st%d%nik) = 0.
!       do j=1,scdm%box_size*2+1
!          do k=1,scdm%box_size*2+1
!             do l=1,scdm%box_size*2+1
!                ip = (j-1)*(2*(scdm%box_size*2+1))**2+(k-1)*(2*(scdm%box_size*2+1)) + l
!                scdm%st%X(psi)(scdm%box(j,k,l,v),st%d%dim,v,scdm%st%d%nik) = scdm%X(psi)(ip,count)
!             enddo
!          enddo
!       enddo
!       !
!    enddo
!!--------------------------------------------------------------------------
!print *, 'orthonrmality in boxes after re-ortho: ================'
!do j=1,nval
!   do i=1,nval
!      print *, i,j, dot_product(scdm%st%X(psi)(1:mesh%np,1,i,1),scdm%st%X(psi)(1:mesh%np,1,j,1))*mesh%volume_element
!   enddo
!enddo
!print *, '======================================================'
!endif
!

! check span
!SAFEx_ALLOCATE(state_global(1:mesh%np_global))
!state_global(:) = M_ZERO
!print *, 'quality of projector:'
!do j=1,nval
!   state_global(:) = M_ZERO
!   do i=1,nval
!      state_global(:) = state_global(:) + &
!           dot_product(scdm%st%X(psi)(:,1,i,1),st%X(psi)(1:mesh%np,1,j,1))*scdm%st%X(psi)(:,1,i,1)*mesh%volume_element
!   enddo
!   print *, j, M_ONE - X(mf_nrm2)(mesh, state_global)
!enddo

    !
    ! set flag to do this only once
    scdm_is_local = .true.
    !
    !
    SAFE_DEALLOCATE_A(Pcc)
    SAFE_DEALLOCATE_A(KSt)
    SAFE_DEALLOCATE_A(scdm_temp)
    !
    !
    call cpu_time(t1)
    if(scdm%root.and.scdm%verbose) call messages_print_var_value(stdout,  'time: all SCDM',t1-t0)
    !
!    return
    !
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
 !            rho(:) = scdm%st%X(psi)(:,st%d%dim,i,scdm%st%d%nik)*scdm%st%X(psi)(:,st%d%dim,j,scdm%st%d%nik)
!!rho(:) = st%X(psi)(:,st%d%dim,i,scdm%st%d%nik)*st%X(psi)(:,st%d%dim,j,scdm%st%d%nik)
 !            rho2(:) = rho(:)
 !            pot(:) = 0.
 !            call X(poisson_solve)(psolver, pot, rho, all_nodes = .false.)
 !            !
 !            ! catch diagonals for double counting
 !            if(i.ne.j) then
 !               exx = exx - 0.5*dot_product(pot(:),rho2(:))*mesh%volume_element
 !            else
 !               exx = exx - 0.25*dot_product(pot(:),rho2(:))*mesh%volume_element
 !            endif
!!          endif
 !         !
 !      enddo
 !   enddo
 !   !
 !   call messages_print_var_value(stdout,'exx[eV] = ', exx*27.211396132)
 !   SxAFE_DEALLOCATE_A(rho)
 !   SxAFE_DEALLOCATE_A(rho2)
 !   SxAFE_DEALLOCATE_A(pot)
    !
    POP_SUB(X(scdm_localize))
    !
  end subroutine X(scdm_localize)


 subroutine X(invert)(n,A)
    ! stupid routine to invert with LAPACK
    ! is very redundant here, shoudl be replaced by something smart
    integer         ::  n
    R_TYPE          ::  A(n,n)
    integer         :: ierror,ipiv(n), lwork
    R_TYPE,pointer  :: work(:)
    FLOAT           :: temp
    !
    PUSH_SUB(X(invert))
    call X(getrf)( n, n, A, n, ipiv, ierror )
    if( ierror.eq.0 ) then
       !workspace query
       call X(getri)( n, A, n, ipiv, temp, -1, ierror )
       lwork = temp ! dimension of workspace
       allocate(work(lwork*2))
       call X(getri)( n, A, n, ipiv, work, lwork, ierror )
    else
       stop 'Terminating due to failed LU decomp'
    end if
    if (ierror.ne.0) then
       stop 'Terminating due to failed inversion'
    end if
    deallocate(work)
    !
    POP_SUB(X(invert))
    !
  end subroutine X(invert)
