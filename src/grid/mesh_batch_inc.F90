!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Verstraete
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

subroutine X(mesh_batch_dotp_matrix)(mesh, aa, bb, dot, symm, reduce)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(in)    :: aa
  type(batch_t),     intent(in)    :: bb
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: symm         !for the moment it is ignored
  logical, optional, intent(in)    :: reduce

  integer :: ist, jst, idim, sp, block_size, ep, ip, ldaa, ldbb, indb, jndb, eff_size
  R_TYPE :: ss, tmp1, tmp2
  R_TYPE, allocatable :: dd(:, :)
#ifdef HAVE_MPI
  R_TYPE, allocatable :: ddtmp(:, :)
#endif
  type(profile_t), save :: prof, profgemm, profcomm
  logical :: use_blas, reduce_, conj
#ifdef HAVE_OPENCL
  type(opencl_mem_t) :: dot_buffer
  type(c_ptr)        :: kernel
#endif

  PUSH_SUB(X(mesh_batch_dotp_matrix))
  call profiling_in(prof, "DOTP_BATCH")

#ifdef HAVE_MPI
  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce
#endif
  conj = .false.

  ASSERT(aa%dim == bb%dim)
  ASSERT(batch_status(aa) == batch_status(bb))

  SAFE_ALLOCATE(dd(1:aa%nst, 1:bb%nst))

  select case(batch_status(aa))
  case(BATCH_NOT_PACKED)
    use_blas = associated(aa%X(psicont)) .and. associated(bb%X(psicont)) .and. (.not. mesh%use_curvilinear) .and. (aa%dim == 1)

    if(use_blas) then
      call profiling_in(profgemm, "DOTP_BATCH_GEMM")

      ldaa = size(aa%X(psicont), dim = 1)
      ldbb = size(bb%X(psicont), dim = 1)
      call blas_gemm('c', 'n', aa%nst, bb%nst, mesh%np, &
        R_TOTYPE(mesh%volume_element), &
        aa%X(psicont)(1, 1, 1), ldaa, &
        bb%X(psicont)(1, 1, 1), ldbb, &
        R_TOTYPE(M_ZERO), dd(1, 1), aa%nst)

    else

      dd = R_TOTYPE(M_ZERO)

      block_size = hardware%X(block_size)

      do idim = 1, aa%dim
        do sp = 1, mesh%np, block_size
          ep = min(mesh%np, sp + block_size - 1)

          if(mesh%use_curvilinear) then

            do ist = 1, aa%nst
              indb = batch_linear_index(aa, (/ist, idim/))
              do jst = 1, bb%nst
                jndb = batch_linear_index(bb, (/jst, idim/))

                ss = M_ZERO
                do ip = sp, ep
                  ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%states_linear(indb)%X(psi)(ip))*bb%states_linear(jndb)%X(psi)(ip)
                end do
                dd(ist, jst) = dd(ist, jst) + ss

              end do
            end do

          else

            do ist = 1, aa%nst
              indb = batch_linear_index(aa, (/ist, idim/))
              do jst = 1, bb%nst
                jndb = batch_linear_index(bb, (/jst, idim/))

                dd(ist, jst) = dd(ist, jst) + mesh%volume_element*&
                  blas_dot(ep - sp + 1, aa%states_linear(indb)%X(psi)(sp), 1, bb%states_linear(jndb)%X(psi)(sp), 1)
              end do
            end do

          end if
        end do
      end do

    end if
  case(BATCH_PACKED)
    ASSERT(.not. mesh%use_curvilinear)
    use_blas = aa%dim == 1

    if(use_blas) then
      conj = .true.
      call profiling_in(profgemm, "DOTP_BATCH_GEMM")

      ldaa = aa%pack%size(1)
      ldbb = bb%pack%size(1)
      call blas_gemm(transa = 'n', transb = 'c', m = aa%nst, n = bb%nst, k = mesh%np, &
        alpha = R_TOTYPE(mesh%volume_element), &
        a = aa%pack%X(psi)(1, 1), lda = ldaa, &
        b = bb%pack%X(psi)(1, 1), ldb = ldbb, &
        beta = R_TOTYPE(M_ZERO), c = dd(1, 1), ldc = aa%nst)
    else

      do ist = 1, aa%nst
        do jst = 1, bb%nst
          tmp1 = 0.0
          tmp2 = 0.0
          do ip = 1, mesh%np
            tmp1 = tmp1 + R_CONJ(aa%pack%X(psi)(2*ist - 1, ip))*bb%pack%X(psi)(2*jst - 1, ip)
            tmp2 = tmp2 + R_CONJ(aa%pack%X(psi)(2*ist    , ip))*bb%pack%X(psi)(2*jst    , ip)
          end do
          dd(ist, jst) = mesh%volume_element*(tmp1 + tmp2)
        end do
      end do

    end if

  case(BATCH_CL_PACKED)
    use_blas = .false.
    ASSERT(.not. mesh%use_curvilinear)
#ifdef HAVE_OPENCL

    kernel = X(kernel_dot_matrix)
#ifdef R_TCOMPLEX
    if(aa%dim > 1) then
      kernel = zkernel_dot_matrix_spinors
    end if
#else
    ASSERT(aa%dim == 1)
#endif

    call opencl_create_buffer(dot_buffer, CL_MEM_WRITE_ONLY, R_TYPE_VAL, aa%nst*bb%nst)
    
    call opencl_set_kernel_arg(kernel, 0, mesh%np)
    call opencl_set_kernel_arg(kernel, 1, aa%pack%buffer)
    call opencl_set_kernel_arg(kernel, 2, log2(aa%pack%size(1)/aa%dim))
    call opencl_set_kernel_arg(kernel, 3, bb%pack%buffer)
    call opencl_set_kernel_arg(kernel, 4, log2(bb%pack%size(1)/aa%dim))
    call opencl_set_kernel_arg(kernel, 5, dot_buffer)
    call opencl_set_kernel_arg(kernel, 6, aa%nst)
    
    call opencl_kernel_run(kernel, (/aa%pack%size(1)/aa%dim, bb%nst/), &
      (/aa%pack%size(1)/aa%dim, 1/))
    
    call opencl_finish()

    call opencl_read_buffer(dot_buffer, aa%nst*bb%nst, dd)
    call opencl_release_buffer(dot_buffer)

    forall(ist = 1:aa%nst, jst = 1:bb%nst) dd(ist, jst) = mesh%volume_element*dd(ist, jst)
#endif

  end select

  if(mesh%use_curvilinear) then
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + 2*R_MUL))
  else
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + R_MUL))
  end if

  if(use_blas) call profiling_out(profgemm)

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "DOTP_BATCH_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dd, dim = (/aa%nst, bb%nst/))
    call profiling_out(profcomm)
  end if
#endif

  if(conj) then
    forall(ist = 1:aa%nst, jst = 1:bb%nst) dot(aa%states(ist)%ist, bb%states(jst)%ist) = R_CONJ(dd(ist, jst))
  else
    forall(ist = 1:aa%nst, jst = 1:bb%nst) dot(aa%states(ist)%ist, bb%states(jst)%ist) = dd(ist, jst)
  end if

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_matrix))
end subroutine X(mesh_batch_dotp_matrix)

!-----------------------------------------------------------------

subroutine X(mesh_batch_dotp_self)(mesh, aa, dot, reduce)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(in)    :: aa
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: reduce

  integer :: ist, jst, idim, sp, block_size, ep, ip, lda, indb, jndb
  R_TYPE :: ss
  type(profile_t), save :: prof, profgemm, profcomm
  logical :: use_blas, reduce_
  R_TYPE, allocatable :: dd(:, :)

  PUSH_SUB(X(mesh_batch_dotp_self))

  ! some limitations of the current implementation
  ASSERT(ubound(dot, dim = 1) >= aa%nst .and. ubound(dot, dim = 2) >= aa%nst)

  if(batch_status(aa) /= BATCH_NOT_PACKED) then
    call X(mesh_batch_dotp_matrix)(mesh, aa, aa, dot, reduce)
    POP_SUB(X(mesh_batch_dotp_self))
    return
  end if

  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce

  use_blas = associated(aa%X(psicont)) .and. (.not. mesh%use_curvilinear)

  SAFE_ALLOCATE(dd(1:aa%nst, 1:aa%nst))

  call profiling_in(prof, "BATCH_DOTP_SELF")

  if(use_blas) then
    call profiling_in(profgemm, "BATCH_HERK")

    ! For some reason this has to be set to zero by hand (a bug in
    ! some Blas libraries?). Otherwise NaNs might contaminate the
    ! result.
    dd(1:aa%nst, 1:aa%nst) = R_TOTYPE(CNST(0.0))

    lda = size(aa%X(psicont), dim = 1)*aa%dim

    call blas_herk('l', 'c', aa%nst, mesh%np, mesh%vol_pp(1), aa%X(psicont)(1, 1, 1), &
      lda, M_ZERO, dd(1, 1), ubound(dd, dim = 1))

    if(aa%dim == 2) then
      call blas_herk('l', 'c', aa%nst, mesh%np, mesh%vol_pp(1), aa%X(psicont)(1, 2, 1), &
        lda, M_ONE, dd(1, 1), ubound(dd, dim = 1))
    end if

  else

    dd = R_TOTYPE(M_ZERO)

    block_size = hardware%X(block_size)

    do idim = 1, aa%dim
      do sp = 1, mesh%np, block_size
        ep = min(mesh%np, sp + block_size - 1)

        if(mesh%use_curvilinear) then

          do ist = 1, aa%nst
            indb = batch_linear_index(aa, (/ist, idim/))
            do jst = 1, ist
              jndb = batch_linear_index(aa, (/jst, idim/))
              ss = M_ZERO
              do ip = sp, ep
                ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%states_linear(indb)%X(psi)(ip))*aa%states_linear(jndb)%X(psi)(ip)
              end do
              dd(ist, jst) = dd(ist, jst) + ss

            end do
          end do

        else

          do ist = 1, aa%nst
            indb = batch_linear_index(aa, (/ist, idim/))
            do jst = 1, ist
              jndb = batch_linear_index(aa, (/jst, idim/))
              dd(ist, jst) = dd(ist, jst) + mesh%volume_element*&
                blas_dot(ep - sp + 1, aa%states_linear(indb)%X(psi)(sp), 1, aa%states_linear(jndb)%X(psi)(sp), 1)
            end do
          end do

        end if
      end do
    end do
  end if

  if(mesh%use_curvilinear) then
    call profiling_count_operations(dble(mesh%np)*aa%nst**2*(R_ADD + 2*R_MUL))
  else
    call profiling_count_operations(dble(mesh%np)*aa%nst**2*(R_ADD + R_MUL))
  end if

  if(use_blas) call profiling_out(profgemm)

  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "BATCH_SELF_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dd, dim = (/aa%nst, aa%nst/))
    call profiling_out(profcomm)
  end if

  forall(ist = 1:aa%nst)
    forall(jst = 1:aa%nst) 
      dot(aa%states(ist)%ist, aa%states(jst)%ist) = dd(ist, jst)
      dot(aa%states(jst)%ist, aa%states(ist)%ist) = R_CONJ(dd(ist, jst))
    end forall
  end forall

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_self))
end subroutine X(mesh_batch_dotp_self)

!-----------------------------------------------------------------

subroutine X(mesh_batch_rotate)(mesh, aa, transf)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(inout) :: aa
  R_TYPE,            intent(in)    :: transf(:, :)

  R_TYPE, allocatable :: psinew(:, :), psicopy(:, :)
  
  integer :: ist, idim, block_size, size, sp, indb
  type(profile_t), save :: prof

  call profiling_in(prof, "BATCH_ROTATE")
  ASSERT(batch_status(aa) == BATCH_NOT_PACKED)

#ifdef R_TREAL  
  block_size = max(40, hardware%l2%size/(2*8*aa%nst))
#else
  block_size = max(20, hardware%l2%size/(2*16*aa%nst))
#endif

  SAFE_ALLOCATE(psinew(1:block_size, 1:aa%nst))
  SAFE_ALLOCATE(psicopy(1:block_size, 1:aa%nst))

  do sp = 1, mesh%np, block_size
    size = min(block_size, mesh%np - sp + 1)
    
    do idim = 1, aa%dim

      do ist = 1, aa%nst
        indb = batch_linear_index(aa, (/ist, idim/))
        call blas_copy(size, aa%states_linear(indb)%X(psi)(sp), 1, psicopy(1, ist), 1)
      end do
      
      call blas_gemm('N', 'N', &
        size, aa%nst, aa%nst, &
        R_TOTYPE(M_ONE), psicopy(1, 1), block_size, &
        transf(1, 1), aa%nst, &
        R_TOTYPE(M_ZERO), psinew(1, 1), block_size)
      
      do ist = 1, aa%nst
        indb = batch_linear_index(aa, (/ist, idim/))
        call blas_copy(size, psinew(1, ist), 1, aa%states_linear(indb)%X(psi)(sp), 1)
      end do
      
    end do
  end do

  SAFE_DEALLOCATE_A(psicopy)
  SAFE_DEALLOCATE_A(psinew)

  call profiling_count_operations((R_ADD + R_MUL)*dble(mesh%np)*aa%dim*aa%nst**2)

  call profiling_out(prof)

end subroutine X(mesh_batch_rotate)

! --------------------------------------------------------------------------

subroutine X(mesh_batch_dotp_vector)(mesh, aa, bb, dot, reduce)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(in)    :: aa
  type(batch_t),     intent(in)    :: bb
  R_TYPE,            intent(inout) :: dot(:)
  logical, optional, intent(in)    :: reduce

  integer :: ist, indb, idim, ip
  logical :: reduce_
  type(profile_t), save :: prof, profcomm
  R_TYPE, allocatable :: tmp(:)
  type(cl_kernel_t), save :: kernel
  type(c_ptr)             :: kernel_ref
#ifdef HAVE_OPENCL
  type(opencl_mem_t)  :: dot_buffer
#endif

  PUSH_SUB(X(mesh_batch_dotp_vector))
  call profiling_in(prof, "DOTPV_BATCH")

  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce
  
  ASSERT(aa%nst == bb%nst)
  ASSERT(aa%dim == bb%dim)

  select case(batch_status(aa))
  case(BATCH_NOT_PACKED)
    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = batch_linear_index(aa, (/ist, idim/))
        dot(ist) = dot(ist) + X(mf_dotp)(mesh, aa%states_linear(indb)%X(psi), bb%states_linear(indb)%X(psi), reduce = .false.)
      end do
    end do

  case(BATCH_PACKED)

    SAFE_ALLOCATE(tmp(1:aa%nst_linear))

    tmp = M_ZERO
    
    if(mesh%use_curvilinear) then
      do ip = 1, mesh%np
        do ist = 1, aa%nst_linear
          tmp(ist) = tmp(ist) + mesh%vol_pp(ip)*R_CONJ(aa%pack%X(psi)(ist, ip))*bb%pack%X(psi)(ist, ip)
        end do
      end do
    else
      do ip = 1, mesh%np
        do ist = 1, aa%nst_linear
          tmp(ist) = tmp(ist) + R_CONJ(aa%pack%X(psi)(ist, ip))*bb%pack%X(psi)(ist, ip)
        end do
      end do
    end if

    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = batch_linear_index(aa, (/ist, idim/))
        dot(ist) = dot(ist) + mesh%volume_element*tmp(indb)
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)

  case(BATCH_CL_PACKED)
    SAFE_ALLOCATE(tmp(1:aa%nst_linear))
#ifdef HAVE_OPENCL
    
    call opencl_create_buffer(dot_buffer, CL_MEM_WRITE_ONLY, R_TYPE_VAL, aa%pack%size(1))

    call cl_kernel_start_call(kernel, 'dot_vector.cl', TOSTRING(X(dot_vector)))

    kernel_ref = cl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, mesh%np)
    call opencl_set_kernel_arg(kernel_ref, 1, aa%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(aa%pack%size(1)))
    call opencl_set_kernel_arg(kernel_ref, 3, bb%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 4, log2(bb%pack%size(1)))
    call opencl_set_kernel_arg(kernel_ref, 5, dot_buffer)
    
    call opencl_kernel_run(kernel_ref, (/aa%pack%size(1)/), (/aa%pack%size(1)/))
    
    call opencl_finish()

    call opencl_read_buffer(dot_buffer, aa%nst_linear, tmp)
    call opencl_release_buffer(dot_buffer)
#endif
    
    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = batch_linear_index(aa, (/ist, idim/))
        dot(ist) = dot(ist) + mesh%volume_element*tmp(indb)
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
  end select

  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "DOTPV_BATCH_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dot, aa%nst)
    call profiling_out(profcomm)
  end if

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_vector))
end subroutine X(mesh_batch_dotp_vector)

!--------------------------------------------------------------------------------------

!> This functions exchanges points of a mesh according to a certain
!! map. Two possible maps can be given. Only one map argument must be present.

subroutine X(mesh_batch_exchange_points)(mesh, aa, forward_map, backward_map)
  type(mesh_t),      intent(in)    :: mesh            !< The mesh descriptor.
  type(batch_t),     intent(inout) :: aa              !< A batch which contains the mesh functions whose points will be exchanged.
  integer, optional, intent(in)    :: forward_map(:)  !< A map which gives the destination of the value each point.
  integer, optional, intent(in)    :: backward_map(:) !< A map which gives the source of the value of each point.

#ifdef HAVE_MPI
  integer :: ip, ipg, npart, ipart, ist, pos, nstl
  integer, allocatable :: send_count(:), recv_count(:), send_disp(:), recv_disp(:)
  R_TYPE, allocatable  :: send_buffer(:, :), recv_buffer(:, :)
#endif

  PUSH_SUB(X(mesh_batch_exchange_points))

  ASSERT(present(backward_map) .neqv. present(forward_map))
  ASSERT(batch_type(aa) == R_TYPE_VAL)
  ASSERT(batch_status(aa) == BATCH_NOT_PACKED)

  if(.not. mesh%parallel_in_domains) then
    message(1) = "Not implemented for the serial case. Really, only in parallel."
    call messages_fatal(1)
  else

#ifdef HAVE_MPI
    npart = mesh%mpi_grp%size
    nstl = aa%nst_linear

    SAFE_ALLOCATE(send_count(1:npart))
    SAFE_ALLOCATE(recv_count(1:npart))
    SAFE_ALLOCATE(send_disp(1:npart))
    SAFE_ALLOCATE(recv_disp(1:npart))
    SAFE_ALLOCATE(send_buffer(1:nstl, mesh%np))

    if(present(forward_map)) then

      ASSERT(ubound(forward_map, dim = 1) == mesh%np_global)

      send_count = 0
      do ip = 1, mesh%np
        !get the global point
        ipg = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
        !the destination
        ipart = mesh%vp%part(forward_map(ipg))
        INCR(send_count(ipart), 1)
      end do

      ASSERT(sum(send_count) == mesh%np)

      recv_count = 0
      do ipg = 1, mesh%np_global
        if(mesh%vp%part(forward_map(ipg)) == mesh%vp%partno) then
          INCR(recv_count(mesh%vp%part(ipg)), 1)
        end if
      end do

      ASSERT(sum(recv_count) == mesh%np)

      send_disp(1) = 0
      recv_disp(1) = 0
      do ipart = 2, npart
        send_disp(ipart) = send_disp(ipart - 1) + send_count(ipart - 1)
        recv_disp(ipart) = recv_disp(ipart - 1) + recv_count(ipart - 1)
      end do

      ASSERT(send_disp(npart) + send_count(npart) == mesh%np)
      ASSERT(recv_disp(npart) + recv_count(npart) == mesh%np)

      !pack for sending
      send_count = 0
      do ip = 1, mesh%np
        !get the global point
        ipg = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
        !the destination
        ipart = mesh%vp%part(forward_map(ipg))
        INCR(send_count(ipart), 1)
        pos = send_disp(ipart) + send_count(ipart)
        forall(ist = 1:nstl) send_buffer(ist, pos) = aa%states_linear(ist)%X(psi)(ip)
      end do

      SAFE_ALLOCATE(recv_buffer(1:nstl, mesh%np))

      call MPI_Alltoallv(send_buffer(1, 1), send_count*nstl, send_disp*nstl, R_MPITYPE, &
        recv_buffer(1, 1), recv_count*nstl, recv_disp*nstl, R_MPITYPE, mesh%mpi_grp%comm, mpi_err)

      SAFE_DEALLOCATE_A(send_buffer)

      recv_count = 0
      do ipg = 1, mesh%np_global
        if(mesh%vp%part(forward_map(ipg)) == mesh%vp%partno) then
          ip = vec_global2local(mesh%vp, forward_map(ipg), mesh%vp%partno)
          ASSERT(ip /= 0)
          ipart = mesh%vp%part(ipg)
          INCR(recv_count(ipart), 1)
          pos = recv_disp(ipart) + recv_count(ipart)
          forall(ist = 1:nstl) aa%states_linear(ist)%X(psi)(ip) = recv_buffer(ist, pos)
        end if
      end do

    else ! backward map

      ASSERT(ubound(backward_map, dim = 1) == mesh%np_global)

      recv_count = 0
      do ip = 1, mesh%np
        !get the global point
        ipg = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
        !the source
        ipart = mesh%vp%part(backward_map(ipg))
        INCR(recv_count(ipart), 1)
      end do

      ASSERT(sum(recv_count) == mesh%np)

      send_count = 0
      do ipg = 1, mesh%np_global
        if(mesh%vp%part(backward_map(ipg)) == mesh%vp%partno) then
          INCR(send_count(mesh%vp%part(ipg)), 1)
        end if
      end do

      ASSERT(sum(send_count) == mesh%np)

      send_disp(1) = 0
      recv_disp(1) = 0
      do ipart = 2, npart
        send_disp(ipart) = send_disp(ipart - 1) + send_count(ipart - 1)
        recv_disp(ipart) = recv_disp(ipart - 1) + recv_count(ipart - 1)
      end do

      ASSERT(send_disp(npart) + send_count(npart) == mesh%np)
      ASSERT(recv_disp(npart) + recv_count(npart) == mesh%np)

      !pack for sending
      send_count = 0
      do ipg = 1, mesh%np_global
        if(mesh%vp%part(backward_map(ipg)) == mesh%vp%partno) then
          ip = vec_global2local(mesh%vp, backward_map(ipg), mesh%vp%partno)
          ipart = mesh%vp%part(ipg)
          INCR(send_count(ipart), 1)
          pos = send_disp(ipart) + send_count(ipart)
          forall(ist = 1:nstl) send_buffer(ist, pos) = aa%states_linear(ist)%X(psi)(ip) 
        end if
      end do

      SAFE_ALLOCATE(recv_buffer(1:nstl, mesh%np))

      call MPI_Alltoallv(send_buffer(1, 1), send_count*nstl, send_disp*nstl, R_MPITYPE, &
        recv_buffer(1, 1), recv_count*nstl, recv_disp*nstl, R_MPITYPE, mesh%mpi_grp%comm, mpi_err)

      SAFE_DEALLOCATE_A(send_buffer)

      recv_count = 0
      do ip = 1, mesh%np
        !get the global point
        ipg = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
        !the destination
        ipart = mesh%vp%part(backward_map(ipg))
        INCR(recv_count(ipart), 1)
        pos = recv_disp(ipart) + recv_count(ipart)
        forall(ist = 1:nstl) aa%states_linear(ist)%X(psi)(ip) = recv_buffer(ist, pos)
      end do

    end if

    SAFE_DEALLOCATE_A(send_count)
    SAFE_DEALLOCATE_A(recv_count)
    SAFE_DEALLOCATE_A(send_disp)
    SAFE_DEALLOCATE_A(recv_disp)
#endif
  end if

  POP_SUB(X(mesh_batch_exchange_points))
end subroutine X(mesh_batch_exchange_points)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
