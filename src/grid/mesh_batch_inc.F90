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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

subroutine X(mesh_batch_dotp_matrix)(mesh, aa, bb, dot, symm, reduce)
  type(mesh_t),      intent(in)    :: mesh
  class(batch_t),    intent(in)    :: aa
  class(batch_t),    intent(in)    :: bb
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: symm         !< for the moment it is ignored
  logical, optional, intent(in)    :: reduce

  integer :: ist, jst, idim, sp, block_size, ep, ip, ldaa, ldbb, indb, jndb
  R_TYPE :: ss, tmp1, tmp2
  R_TYPE, allocatable :: dd(:, :)
  logical :: use_blas, conj
  type(accel_mem_t) :: dot_buffer
  type(profile_t), save :: prof_copy, prof_gemmcl, prof, profgemm
  integer :: wgsize
  integer :: local_sizes(3)
  integer :: global_sizes(3)

  logical :: reduce_
  type(profile_t), save :: profcomm
  
  PUSH_SUB(X(mesh_batch_dotp_matrix))
  call profiling_in(prof, "DOTP_BATCH")

  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce
  conj = .false.

  call aa%check_compatibility_with(bb, only_check_dim = .true.)

  SAFE_ALLOCATE(dd(1:aa%nst, 1:bb%nst))
  ! This has to be set to zero by hand since NaN * 0 = NaN.
  dd(1:aa%nst, 1:bb%nst) = R_TOTYPE(M_ZERO)

  use_blas = .false.
  
  select case(aa%status())
  case(BATCH_NOT_PACKED)
    use_blas = associated(aa%X(ff)) .and. associated(bb%X(ff)) .and. (.not. mesh%use_curvilinear) .and. (aa%dim == 1)

    if(use_blas) then
      call profiling_in(profgemm, "DOTP_BATCH_GEMM")

      ldaa = size(aa%X(ff), dim = 1)
      ldbb = size(bb%X(ff), dim = 1)

      call lalg_gemmt(aa%nst, aa%dim, bb%nst, bb%dim, mesh%np, R_TOTYPE(mesh%volume_element), &
        aa%X(ff), bb%X(ff), R_TOTYPE(M_ZERO), dd)

    else

      block_size = hardware%X(block_size)

      do idim = 1, aa%dim
        do sp = 1, mesh%np, block_size
          ep = min(mesh%np, sp + block_size - 1)

          if(mesh%use_curvilinear) then

            do ist = 1, aa%nst
              indb = aa%ist_idim_to_linear((/ist, idim/))
              do jst = 1, bb%nst
                jndb = bb%ist_idim_to_linear((/jst, idim/))

                ss = M_ZERO
                do ip = sp, ep
                  ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%X(ff_linear)(ip, indb))*bb%X(ff_linear)(ip, jndb)
                end do
                dd(ist, jst) = dd(ist, jst) + ss

              end do
            end do

          else

            do ist = 1, aa%nst
              indb = aa%ist_idim_to_linear((/ist, idim/))
              do jst = 1, bb%nst
                jndb = bb%ist_idim_to_linear((/jst, idim/))

                dd(ist, jst) = dd(ist, jst) + mesh%volume_element*&
                  blas_dot(ep - sp + 1, aa%X(ff_linear)(sp, indb), 1, bb%X(ff_linear)(sp, jndb), 1)
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

      ldaa = aa%pack_size(1)
      ldbb = bb%pack_size(1)
      call blas_gemm(transa = 'n', transb = 'c', m = aa%nst, n = bb%nst, k = mesh%np, &
        alpha = R_TOTYPE(mesh%volume_element), &
        a = aa%X(ff_pack)(1, 1), lda = ldaa, &
        b = bb%X(ff_pack)(1, 1), ldb = ldbb, &
        beta = R_TOTYPE(M_ZERO), c = dd(1, 1), ldc = aa%nst)
      
    else

      do ist = 1, aa%nst
        do jst = 1, bb%nst
          tmp1 = M_ZERO
          tmp2 = M_ZERO
          do ip = 1, mesh%np
            tmp1 = tmp1 + R_CONJ(aa%X(ff_pack)(2*ist - 1, ip))*bb%X(ff_pack)(2*jst - 1, ip)
            tmp2 = tmp2 + R_CONJ(aa%X(ff_pack)(2*ist    , ip))*bb%X(ff_pack)(2*jst    , ip)
          end do
          dd(ist, jst) = mesh%volume_element*(tmp1 + tmp2)
        end do
      end do

    end if

  case(BATCH_DEVICE_PACKED)
    ASSERT(.not. mesh%use_curvilinear)

    if(aa%dim==1) then
 
      call accel_create_buffer(dot_buffer, ACCEL_MEM_WRITE_ONLY, R_TYPE_VAL, aa%nst*bb%nst)

      call profiling_in(prof_gemmcl, "DOTP_BATCH_CL_GEMM")
      
      call X(accel_gemm)(transA = CUBLAS_OP_N, transB = CUBLAS_OP_T, &
        M = int(aa%nst, 8), N = int(bb%nst, 8), K = int(mesh%np, 8), alpha = R_TOTYPE(M_ONE), &
        A = aa%ff_device, offA = 0_8, lda = int(aa%pack_size(1), 8), &
        B = bb%ff_device, offB = 0_8, ldb = int(bb%pack_size(1), 8), beta = R_TOTYPE(M_ZERO), &
        C = dot_buffer, offC = 0_8, ldc = int(aa%nst, 8))
  
      call profiling_count_operations(TOFLOAT(mesh%np)*aa%nst*bb%nst*(R_ADD + R_MUL))
  
      call accel_finish()
      call profiling_out(prof_gemmcl)
  
      call profiling_in(prof_copy, 'DOTP_BATCH_COPY')
      call accel_read_buffer(dot_buffer, aa%nst*bb%nst, dd)
      call profiling_count_transfers(aa%nst*bb%nst, dd(1, 1))
      call accel_finish()
      call profiling_out(prof_copy)
  
      call accel_release_buffer(dot_buffer)

    else

      ASSERT(R_TYPE_VAL == TYPE_CMPLX)

      call accel_create_buffer(dot_buffer, ACCEL_MEM_WRITE_ONLY, R_TYPE_VAL, aa%nst*bb%nst)

      wgsize = accel_kernel_workgroup_size(zkernel_dot_matrix_spinors)

      global_sizes = (/ pad(aa%nst, wgsize/bb%nst),  bb%nst, 1 /)
      local_sizes  = (/ wgsize/bb%nst,               bb%nst, 1 /)
     
      ASSERT(accel_buffer_is_allocated(aa%ff_device))
      ASSERT(accel_buffer_is_allocated(bb%ff_device))
      ASSERT(accel_buffer_is_allocated(dot_buffer))

      call profiling_in(prof_gemmcl, "DOTP_BATCH_CL_KERNEL")

      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 0, mesh%np)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 1, aa%nst)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 2, bb%nst)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 3, aa%ff_device)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 4, log2(aa%pack_size(1)))
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 5, bb%ff_device)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 6, log2(bb%pack_size(1)))
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 7, dot_buffer)
      call accel_set_kernel_arg(zkernel_dot_matrix_spinors, 8, aa%nst)


      call accel_kernel_run(zkernel_dot_matrix_spinors, global_sizes, local_sizes)
  
      call accel_finish()
      call profiling_count_operations(TOFLOAT(aa%nst*bb%nst*(mesh%np*(R_ADD + R_MUL)) + R_ADD )) ! check !!


      call profiling_out(prof_gemmcl)
  
      call profiling_in(prof_copy, 'DOTP_BATCH_COPY')
      call accel_read_buffer(dot_buffer, aa%nst*bb%nst, dd)
      call profiling_count_transfers(aa%nst*bb%nst, dd(1, 1))
      call accel_finish()
      call profiling_out(prof_copy)
  
      call accel_release_buffer(dot_buffer)

    end if

    do ist = 1, aa%nst
      do jst = 1, bb%nst
        dd(ist, jst) = mesh%volume_element*dd(ist, jst)
      end do
    end do

  case default
    ASSERT(.false.)

  end select

  if(aa%status() /= BATCH_DEVICE_PACKED) then
    if(mesh%use_curvilinear) then
      call profiling_count_operations(TOFLOAT(mesh%np)*aa%nst*bb%nst*aa%dim*(R_ADD + 2*R_MUL))
    else
      call profiling_count_operations(TOFLOAT(mesh%np)*aa%nst*bb%nst*aa%dim*(R_ADD + R_MUL))
    end if
  end if

  if(use_blas) call profiling_out(profgemm)

  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "DOTP_BATCH_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dd)
    call profiling_out(profcomm)
  end if

  if(conj) then
    do jst = 1, bb%nst
      do ist = 1, aa%nst
        dot(aa%ist(ist), bb%ist(jst)) = R_CONJ(dd(ist, jst))
      end do
    end do
  else
    do jst = 1, bb%nst
      do ist = 1, aa%nst
        dot(aa%ist(ist), bb%ist(jst)) = dd(ist, jst)
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_matrix))
end subroutine X(mesh_batch_dotp_matrix)

!-----------------------------------------------------------------

subroutine X(mesh_batch_dotp_self)(mesh, aa, dot, reduce)
  type(mesh_t),      intent(in)    :: mesh
  class(batch_t),    intent(in)    :: aa
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

  if(aa%status() /= BATCH_NOT_PACKED) then
    call X(mesh_batch_dotp_matrix)(mesh, aa, aa, dot, reduce)
    POP_SUB(X(mesh_batch_dotp_self))
    return
  end if

  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce

  use_blas = associated(aa%X(ff)) .and. (.not. mesh%use_curvilinear)

  SAFE_ALLOCATE(dd(1:aa%nst, 1:aa%nst))
  ! This has to be set to zero by hand since NaN * 0 = NaN.
  dd(1:aa%nst, 1:aa%nst) = R_TOTYPE(CNST(0.0))

  call profiling_in(prof, "BATCH_DOTP_SELF")

  if(use_blas) then
    call profiling_in(profgemm, "BATCH_HERK")

    lda = size(aa%X(ff), dim = 1)*aa%dim

    call blas_herk('l', 'c', aa%nst, mesh%np, mesh%vol_pp(1), aa%X(ff)(1, 1, 1), &
      lda, M_ZERO, dd(1, 1), ubound(dd, dim = 1))

    if(aa%dim == 2) then
      call blas_herk('l', 'c', aa%nst, mesh%np, mesh%vol_pp(1), aa%X(ff)(1, 2, 1), &
        lda, M_ONE, dd(1, 1), ubound(dd, dim = 1))
    end if

  else

    block_size = hardware%X(block_size)

    do idim = 1, aa%dim
      do sp = 1, mesh%np, block_size
        ep = min(mesh%np, sp + block_size - 1)

        if(mesh%use_curvilinear) then

          do ist = 1, aa%nst
            indb = aa%ist_idim_to_linear((/ist, idim/))
            do jst = 1, ist
              jndb = aa%ist_idim_to_linear((/jst, idim/))
              ss = M_ZERO
              do ip = sp, ep
                ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%X(ff_linear)(ip, indb))*aa%X(ff_linear)(ip, jndb)
              end do
              dd(ist, jst) = dd(ist, jst) + ss

            end do
          end do

        else

          do ist = 1, aa%nst
            indb = aa%ist_idim_to_linear((/ist, idim/))
            do jst = 1, ist
              jndb = aa%ist_idim_to_linear((/jst, idim/))
              dd(ist, jst) = dd(ist, jst) + mesh%volume_element*&
                blas_dot(ep - sp + 1, aa%X(ff_linear)(sp, indb), 1, aa%X(ff_linear)(sp, jndb), 1)
            end do
          end do

        end if
      end do
    end do
  end if

  if(mesh%use_curvilinear) then
    call profiling_count_operations(TOFLOAT(mesh%np)*aa%nst**2*aa%dim*(R_ADD + 2*R_MUL))
  else
    call profiling_count_operations(TOFLOAT(mesh%np)*aa%nst**2*aa%dim*(R_ADD + R_MUL))
  end if

  if(use_blas) call profiling_out(profgemm)

  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "BATCH_SELF_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dd)
    call profiling_out(profcomm)
  end if

  do ist = 1, aa%nst
    do jst = 1, ist
      dot(aa%ist(ist), aa%ist(jst)) = dd(ist, jst)
      dot(aa%ist(jst), aa%ist(ist)) = R_CONJ(dd(ist, jst))
    end do
  end do

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_self))
end subroutine X(mesh_batch_dotp_self)

! --------------------------------------------------------------------------

subroutine X(mesh_batch_dotp_vector)(mesh, aa, bb, dot, reduce, cproduct)
  type(mesh_t),      intent(in)    :: mesh
  class(batch_t),    intent(in)    :: aa
  class(batch_t),    intent(in)    :: bb
  R_TYPE,            intent(inout) :: dot(:)
  logical, optional, intent(in)    :: reduce
  logical, optional, intent(in)    :: cproduct

  integer :: ist, indb, idim, ip, status
  logical :: cproduct_
  type(profile_t), save :: prof, profcomm
  R_TYPE, allocatable :: tmp(:), cltmp(:, :)
  type(accel_mem_t)  :: dot_buffer

  PUSH_SUB(X(mesh_batch_dotp_vector))
  call profiling_in(prof, "DOTPV_BATCH")

  cproduct_ = optional_default(cproduct, .false.)
  
  call aa%check_compatibility_with(bb)

  status = aa%status()
  ASSERT(bb%status() == status)

  select case(status)
  case(BATCH_NOT_PACKED)
    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        dot(ist) = dot(ist) + X(mf_dotp)(mesh, aa%X(ff_linear)(:, indb), bb%X(ff_linear)(:, indb),& 
           reduce = .false., dotu = cproduct_)
      end do
    end do

  case(BATCH_PACKED)

    SAFE_ALLOCATE(tmp(1:aa%nst_linear))

    tmp = M_ZERO
    
    if(mesh%use_curvilinear) then
      if(.not. cproduct_) then
        !$omp parallel do private(ip, ist) reduction(+:tmp)
        do ip = 1, mesh%np
          do ist = 1, aa%nst_linear
            tmp(ist) = tmp(ist) + mesh%vol_pp(ip)*R_CONJ(aa%X(ff_pack)(ist, ip))*bb%X(ff_pack)(ist, ip)
          end do
        end do
      else
        !$omp parallel do private(ip, ist) reduction(+:tmp)
        do ip = 1, mesh%np
          do ist = 1, aa%nst_linear
            tmp(ist) = tmp(ist) + mesh%vol_pp(ip)*aa%X(ff_pack)(ist, ip)*bb%X(ff_pack)(ist, ip)
          end do
        end do
      end if
    else
      if(.not. cproduct_) then
        !$omp parallel do private(ip, ist) reduction(+:tmp)
        do ip = 1, mesh%np
          do ist = 1, aa%nst_linear
            tmp(ist) = tmp(ist) + R_CONJ(aa%X(ff_pack)(ist, ip))*bb%X(ff_pack)(ist, ip)
          end do
        end do
      else
        !$omp parallel do private(ip, ist) reduction(+:tmp)
        do ip = 1, mesh%np
          do ist = 1, aa%nst_linear
            tmp(ist) = tmp(ist) + aa%X(ff_pack)(ist, ip)*bb%X(ff_pack)(ist, ip)
          end do
        end do
      end if
    end if

    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        dot(ist) = dot(ist) + mesh%volume_element*tmp(indb)
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)

  case(BATCH_DEVICE_PACKED)

    call accel_create_buffer(dot_buffer, ACCEL_MEM_WRITE_ONLY, R_TYPE_VAL, aa%pack_size(1))

    do ist = 1, aa%nst_linear
      call accel_set_stream(ist)
      call X(accel_dot)(n = int(mesh%np, 8), &
        x = aa%ff_device, offx = int(ist - 1, 8), incx = int(aa%pack_size(1), 8), &
        y = bb%ff_device, offy = int(ist - 1, 8), incy = int(bb%pack_size(1), 8), &
        res = dot_buffer, offres = int(ist - 1, 8))
    end do
    call accel_synchronize_all_streams()
    call accel_set_stream(1)

    SAFE_ALLOCATE(cltmp(1:aa%pack_size(1), 1))

    call accel_read_buffer(dot_buffer, aa%pack_size(1), cltmp)

    call accel_release_buffer(dot_buffer)


    do ist = 1, aa%nst
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        dot(ist) = dot(ist) + mesh%volume_element*cltmp(indb, 1)
      end do
    end do

  end select

  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(profcomm, "DOTPV_BATCH_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dot, dim = aa%nst)
    call profiling_out(profcomm)
  end if
  
  call profiling_count_operations(aa%nst_linear*TOFLOAT(mesh%np)*(R_ADD + R_MUL))

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_dotp_vector))
end subroutine X(mesh_batch_dotp_vector)

! --------------------------------------------------------------------------

subroutine X(mesh_batch_mf_dotp)(mesh, aa, psi, dot, reduce, nst)
  type(mesh_t),      intent(in)    :: mesh
  class(batch_t),    intent(in)    :: aa
  R_TYPE,            intent(in)    :: psi(:,:) 
  R_TYPE,            intent(inout) :: dot(:)
  logical, optional, intent(in)    :: reduce
  integer, optional, intent(in)    :: nst

  integer :: ist, indb, idim, ip, nst_
  type(profile_t), save :: prof, profcomm
  R_TYPE, allocatable :: phi(:, :)

  ! Variables related to the GPU:
  type(accel_mem_t) :: psi_buffer
  type(accel_mem_t) :: dot_buffer
  integer :: wgsize, np_padded
  integer :: local_sizes(3)
  integer :: global_sizes(3)

  PUSH_SUB(X(mesh_batch_mf_dotp))
  call profiling_in(prof, "DOTPV_MF_BATCH")

  ASSERT(aa%dim == ubound(psi,dim=2))

  nst_ = aa%nst
  if(present(nst)) nst_ = nst 

  select case(aa%status())
  case(BATCH_NOT_PACKED)
    do ist = 1, nst_
      dot(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        dot(ist) = dot(ist) + X(mf_dotp)(mesh, aa%X(ff_linear)(:, indb), psi(1:mesh%np,idim),& 
           reduce = .false.)
      end do
    end do

  case(BATCH_PACKED)

    SAFE_ALLOCATE(phi(1:mesh%np, aa%dim))

    if(aa%dim == 1) then
      !Here we compute the complex conjuguate of the dot product first and then
      !we take the conjugate at the end

      ! Note: this is to avoid taking the complex conjugate of the whole batch, but rather that of
      ! the single function only.
      ! In the aa%dim>1 case, that is taken care of by the mf_dotp function.

      if(mesh%use_curvilinear) then
        !$omp parallel do
        do ip = 1, mesh%np
          phi(ip, 1) = mesh%vol_pp(ip)*R_CONJ(psi(ip, 1))
        end do
      else
        !$omp parallel do
        do ip = 1, mesh%np
          phi(ip, 1) = R_CONJ(psi(ip, 1))
        end do
      end if

      call blas_gemv('N', nst_, mesh%np, R_TOTYPE(mesh%volume_element), aa%X(ff_pack)(1,1), & 
               ubound(aa%X(ff_pack), dim=1), phi(1,1), 1, R_TOTYPE(M_ZERO), dot(1), 1)

      do ist = 1, nst_
        dot(ist) = R_CONJ(dot(ist))
      end do

    else

      ! Note: curvilinear coordinates are handled inside the mf_dotp function!
  
      dot(1:nst_) = M_ZERO
      do ist = 1, nst_
        call batch_get_state(aa, ist, mesh%np, phi)
        dot(ist) = X(mf_dotp)(mesh, aa%dim, phi(1:mesh%np, 1:aa%dim), psi(1:mesh%np, 1:aa%dim),&
               reduce = .false.)
      end do

    end if

    SAFE_DEALLOCATE_A(phi)

  case(BATCH_DEVICE_PACKED)

    ASSERT(.not. mesh%use_curvilinear)

    np_padded = pad_pow2(mesh%np)

    call accel_create_buffer(dot_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, aa%nst)
    call accel_create_buffer(psi_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, np_padded * aa%dim)

    do idim= 1, aa%dim
      call accel_write_buffer(psi_buffer, mesh%np, psi(1:mesh%np,idim), offset=(idim-1)*np_padded)
    end do
       
    wgsize = accel_kernel_workgroup_size(X(kernel_batch_dotp))

    global_sizes = (/ pad(aa%nst, wgsize),  1, 1 /)
    local_sizes  = (/ wgsize,               1, 1 /)
   
    ASSERT(accel_buffer_is_allocated(aa%ff_device))
    ASSERT(accel_buffer_is_allocated(psi_buffer))
    ASSERT(accel_buffer_is_allocated(dot_buffer))

    call accel_set_kernel_arg(X(kernel_batch_dotp), 0, mesh%np)
    call accel_set_kernel_arg(X(kernel_batch_dotp), 1, nst_)
    call accel_set_kernel_arg(X(kernel_batch_dotp), 2, aa%dim)
    call accel_set_kernel_arg(X(kernel_batch_dotp), 3, aa%ff_device)
    call accel_set_kernel_arg(X(kernel_batch_dotp), 4, log2(aa%pack_size(1)))
    call accel_set_kernel_arg(X(kernel_batch_dotp), 5, psi_buffer)
    call accel_set_kernel_arg(X(kernel_batch_dotp), 6, log2(np_padded))
    call accel_set_kernel_arg(X(kernel_batch_dotp), 7, dot_buffer)

    call accel_kernel_run(X(kernel_batch_dotp), global_sizes, local_sizes)
    call accel_finish() 

    call accel_read_buffer(dot_buffer, nst_, dot)

    call accel_release_buffer(psi_buffer)
    call accel_release_buffer(dot_buffer)

    do ist = 1, nst_
      dot(ist) = dot(ist) * mesh%volume_element
    end do

  end select

  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(profcomm, "DOTPV_MF_BATCH_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dot, dim = nst_)
    call profiling_out(profcomm)
  end if
  
  call profiling_count_operations(nst_*aa%dim*TOFLOAT(mesh%np)*(R_ADD + R_MUL))

  call profiling_out(prof)
  POP_SUB(X(mesh_batch_mf_dotp))
end subroutine X(mesh_batch_mf_dotp)



!--------------------------------------------------------------------------------------

!> This functions exchanges points of a mesh according to a certain
!! map. Two possible maps can be given. Only one map argument must be present.

subroutine X(mesh_batch_exchange_points)(mesh, aa, forward_map, backward_map)
  type(mesh_t),      intent(in)    :: mesh            !< The mesh descriptor.
  class(batch_t),    intent(inout) :: aa              !< A batch which contains the mesh functions whose points will be exchanged.
  integer, optional, intent(in)    :: forward_map(:)  !< A map which gives the destination of the value each point.
  logical, optional, intent(in)    :: backward_map    !< A map which gives the source of the value of each point.
  logical :: packed_on_entry

#ifdef HAVE_MPI
  integer :: ip, ipg, npart, ipart, ist, pos, nstl, np_points, np_inner, np_bndry
  integer, allocatable :: send_count(:), recv_count(:), send_disp(:), recv_disp(:), &
       points_inner(:), points_bndry(:), partno_inner(:), partno_bndry(:)
  integer, allocatable :: send_count_nstl(:), recv_count_nstl(:), send_disp_nstl(:), recv_disp_nstl(:)
  R_TYPE, allocatable  :: send_buffer(:, :), recv_buffer(:, :)
#endif

  PUSH_SUB(X(mesh_batch_exchange_points))

  ASSERT(present(backward_map) .neqv. present(forward_map))
  ASSERT(aa%type() == R_TYPE_VAL)
  packed_on_entry = aa%status() == BATCH_NOT_PACKED
  if (packed_on_entry) then
    call aa%do_unpack(force=.true.)
  end if

  if(.not. mesh%parallel_in_domains) then
    message(1) = "Not implemented for the serial case. Really, only in parallel."
    call messages_fatal(1)
  else

#ifdef HAVE_MPI
    npart = mesh%mpi_grp%size
    nstl = aa%nst_linear

    SAFE_ALLOCATE(send_count(1:npart))
    SAFE_ALLOCATE(recv_count(1:npart))
    SAFE_ALLOCATE(send_count_nstl(1:npart))
    SAFE_ALLOCATE(recv_count_nstl(1:npart))
    SAFE_ALLOCATE(send_disp_nstl(1:npart))
    SAFE_ALLOCATE(recv_disp_nstl(1:npart))
    SAFE_ALLOCATE(send_buffer(1:nstl, mesh%np))
    SAFE_ALLOCATE(recv_buffer(1:nstl, mesh%np))

    if(present(forward_map)) then

      SAFE_ALLOCATE(send_disp(1:npart))
      SAFE_ALLOCATE(recv_disp(1:npart))
      SAFE_ALLOCATE(points_inner(1:mesh%np))
      SAFE_ALLOCATE(points_bndry(1:mesh%np))
      ASSERT(ubound(forward_map, dim = 1) == mesh%np_global)

      send_count = 0
      np_inner   = 0
      np_bndry   = 0
      np_points  = 0
      do ip = 1, mesh%np
        ! Get the temporally global point
        ipg = mesh%vp%local(mesh%vp%xlocal + ip - 1)
        ! Store the global point
        ! Global index can be either in the mesh or in the boundary.
        ! Different treatment is needed for each case.
        if (ipg > mesh%np_global) then
          np_bndry = np_bndry + 1
          points_bndry(np_bndry) = forward_map(ipg) - mesh%np_global
        else
          np_inner = np_inner + 1
          points_inner(np_inner) = forward_map(ipg)
        end if
        np_points = np_points + 1
      end do

      SAFE_ALLOCATE(partno_inner(1:np_points))
      SAFE_ALLOCATE(partno_bndry(1:np_points))
      call partition_get_partition_number(mesh%inner_partition, np_inner, &
           points_inner, partno_inner)
      call partition_get_partition_number(mesh%bndry_partition, np_bndry, &
           points_bndry, partno_bndry)
      SAFE_DEALLOCATE_A(points_inner)
      SAFE_DEALLOCATE_A(points_bndry)
      do ip = 1, np_inner
        ! the destination
        ipart = partno_inner(ip)
        INCR(send_count(ipart), 1)
      end do
      do ip = 1, np_bndry
        ! the destination
        ipart = partno_bndry(ip)
        INCR(send_count(ipart), 1)
      end do
      ASSERT(sum(send_count) == mesh%np)

      ! Receiving number of points is the inverse matrix of the sending points
      call mpi_debug_in(mesh%mpi_grp%comm, C_MPI_ALLTOALL)
      call MPI_Alltoall(send_count(1), 1, MPI_INTEGER, &
                        recv_count(1), 1, MPI_INTEGER, &
                        mesh%mpi_grp%comm, mpi_err)
      call mpi_debug_out(mesh%mpi_grp%comm, C_MPI_ALLTOALL)
      ASSERT(sum(recv_count) == mesh%np)

      send_disp(1) = 0
      recv_disp(1) = 0
      do ipart = 2, npart
        send_disp(ipart) = send_disp(ipart - 1) + send_count(ipart - 1)
        recv_disp(ipart) = recv_disp(ipart - 1) + recv_count(ipart - 1)
      end do

      ASSERT(send_disp(npart) + send_count(npart) == mesh%np)
      ASSERT(recv_disp(npart) + recv_count(npart) == mesh%np)

      ! Pack for sending
      send_count = 0
      ! First inner points
      do ip = 1, np_inner
        !the destination
        ipart = partno_inner(ip)
        INCR(send_count(ipart), 1)
        pos = send_disp(ipart) + send_count(ipart)
        do ist = 1, nstl
          send_buffer(ist, pos) = aa%X(ff_linear)(ip, ist)
        end do
      end do
      ! Then boundary points
      do ip = 1, np_bndry
        !the destination
        ipart = partno_bndry(ip)
        INCR(send_count(ipart), 1)
        pos = send_disp(ipart) + send_count(ipart)
        do ist = 1, nstl
          send_buffer(ist, pos) = aa%X(ff_linear)(ip, ist)
        end do
      end do

      SAFE_DEALLOCATE_A(partno_bndry)
      SAFE_DEALLOCATE_A(partno_inner)

      send_count_nstl = send_count * nstl
      send_disp_nstl = send_disp * nstl
      recv_count_nstl = recv_count * nstl
      recv_disp_nstl = recv_disp * nstl
      call mpi_debug_in(mesh%mpi_grp%comm, C_MPI_ALLTOALLV)
      call MPI_Alltoallv(send_buffer(1, 1), send_count_nstl, send_disp_nstl, R_MPITYPE, &
        recv_buffer(1, 1), recv_count_nstl, recv_disp_nstl, R_MPITYPE, mesh%mpi_grp%comm, mpi_err)
      call mpi_debug_out(mesh%mpi_grp%comm, C_MPI_ALLTOALLV)

      recv_count = 0
      do ipg = 1, mesh%np_global
        if(mesh%vp%part_vec(forward_map(ipg)) == mesh%vp%partno) then
          ip = vec_global2local(mesh%vp, forward_map(ipg), mesh%vp%partno)
          ASSERT(ip /= 0)
          ipart = mesh%vp%part_vec(ipg)
          INCR(recv_count(ipart), 1)
          pos = recv_disp(ipart) + recv_count(ipart)
          do ist = 1, nstl
            aa%X(ff_linear)(ip, ist) = recv_buffer(ist, pos)
          end do
        end if
      end do

      SAFE_DEALLOCATE_A(send_disp)
      SAFE_DEALLOCATE_A(recv_disp)

    else ! backward map

      recv_count = mesh%vp%recv_count
      ASSERT(sum(recv_count) == mesh%np)

      send_count = mesh%vp%send_count
      ASSERT(sum(send_count) == mesh%np)

      ASSERT(mesh%vp%send_disp(npart) + send_count(npart) == mesh%np)
      ASSERT(mesh%vp%recv_disp(npart) + recv_count(npart) == mesh%np)

      ! Pack for sending
      send_count = 0  
      do ip = 1, mesh%np
        ipart = mesh%vp%part_local(ip)
        INCR(send_count(ipart), 1)
        pos = mesh%vp%send_disp(ipart) + send_count(ipart)
        do ist = 1, nstl
          send_buffer(ist, pos) = aa%X(ff_linear)(ip, ist)
        end do
      end do

      send_count_nstl = send_count * nstl
      send_disp_nstl = mesh%vp%send_disp * nstl
      recv_count_nstl = recv_count * nstl
      recv_disp_nstl = mesh%vp%recv_disp * nstl
      call mpi_debug_in(mesh%mpi_grp%comm, C_MPI_ALLTOALLV)
      call MPI_Alltoallv(send_buffer(1, 1), send_count_nstl, send_disp_nstl, R_MPITYPE, &
        recv_buffer(1, 1), recv_count_nstl, recv_disp_nstl, R_MPITYPE, mesh%mpi_grp%comm, mpi_err)
      call mpi_debug_out(mesh%mpi_grp%comm, C_MPI_ALLTOALLV)

      ! Unpack on receiving
      recv_count = 0
      do ip = 1, mesh%np
        ! get the destination
        ipart = mesh%vp%part_local_rev(ip)
        INCR(recv_count(ipart), 1)
        pos = mesh%vp%recv_disp(ipart) + recv_count(ipart)
        do ist = 1, nstl
          aa%X(ff_linear)(ip, ist) = recv_buffer(ist, pos)
        end do
      end do

    end if

    SAFE_DEALLOCATE_A(send_count)
    SAFE_DEALLOCATE_A(recv_count)
    SAFE_DEALLOCATE_A(send_buffer)
    SAFE_DEALLOCATE_A(recv_buffer)
    SAFE_DEALLOCATE_A(send_count_nstl)
    SAFE_DEALLOCATE_A(recv_count_nstl)
    SAFE_DEALLOCATE_A(send_disp_nstl)
    SAFE_DEALLOCATE_A(recv_disp_nstl)
#endif
  end if

  if (packed_on_entry) then
    call aa%do_pack()
  end if
  POP_SUB(X(mesh_batch_exchange_points))
end subroutine X(mesh_batch_exchange_points)

! -----------------------------------------------------
!> This function should not be called directly, but through mesh_batch_nrm2.
subroutine X(priv_mesh_batch_nrm2)(mesh, aa, nrm2)
  type(mesh_t),            intent(in)    :: mesh
  class(batch_t),          intent(in)    :: aa
  FLOAT,                   intent(out)   :: nrm2(:)

  integer :: ist, idim, indb, ip, sp, np, num_threads, ithread
  FLOAT :: a0
  FLOAT, allocatable :: scal(:,:), ssq(:,:)
  type(accel_mem_t)  :: nrm2_buffer
  type(profile_t), save :: prof

  PUSH_SUB(X(priv_mesh_batch_nrm2))
  call profiling_in(prof, 'MESH_BATCH_NRM2')

  select case(aa%status())
  case(BATCH_NOT_PACKED)
    do ist = 1, aa%nst
      nrm2(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        nrm2(ist) = hypot(nrm2(ist), X(mf_nrm2)(mesh, aa%X(ff_linear)(:, indb), reduce = .false.))
      end do
    end do

  case(BATCH_PACKED)
    

    num_threads = 1
    !$omp parallel shared(num_threads)
    !$ num_threads = omp_get_num_threads()
    !$omp end parallel

    SAFE_ALLOCATE(scal(1:aa%nst_linear, 1:num_threads))
    SAFE_ALLOCATE(ssq(1:aa%nst_linear, 1:num_threads))

    scal = M_ZERO
    ssq  = M_ONE

    ! divide the range from 1:mesh%np across the OpenMP threads and sum independently
    ! the reduction is done outside the parallel region
    !$omp parallel private(ithread, sp, np, a0, ip, ist) shared(ssq, scal, num_threads)
    call multicomm_divide_range_omp(mesh%np, sp, np)
    ithread = 1
    !$ ithread = omp_get_thread_num() + 1
    
    if(.not. mesh%use_curvilinear) then

      do ip = sp, sp + np - 1
        do ist = 1, aa%nst_linear
          a0 = abs(aa%X(ff_pack)(ist, ip))
          if(a0 <= M_EPSILON) cycle
          if(scal(ist, ithread) < a0) then
            ssq(ist, ithread) = M_ONE + ssq(ist, ithread)*(scal(ist, ithread)/a0)**2
            scal(ist, ithread) = a0
          else
            ssq(ist, ithread) = ssq(ist, ithread) + (a0/scal(ist, ithread))**2
          end if
        end do
      end do

    else

      do ip = sp, sp + np - 1
        do ist = 1, aa%nst_linear
          a0 = abs(aa%X(ff_pack)(ist, ip))
          if(a0 < M_EPSILON) cycle
          if(scal(ist, ithread) < a0) then
            ssq(ist, ithread) =  mesh%vol_pp(ip)*M_ONE + ssq(ist, ithread)*(scal(ist, ithread)/a0)**2
            scal(ist, ithread) = a0
          else
            ssq(ist, ithread) = ssq(ist, ithread) + mesh%vol_pp(ip)*(a0/scal(ist, ithread))**2
          end if
        end do
      end do

    end if
    !$omp end parallel

    ! now do the reduction: sum the components of the different threads without overflow
    do ithread = 2, num_threads
      do ist = 1, aa%nst_linear
        if (scal(ist, ithread) < M_EPSILON) cycle
        if (scal(ist, 1) < scal(ist, ithread)) then
          ssq(ist, 1) = ssq(ist, 1) * (scal(ist, 1)/scal(ist, ithread))**2 + ssq(ist, ithread)
          scal(ist, 1) = scal(ist, ithread)
        else
          ssq(ist, 1) = ssq(ist, 1) + ssq(ist, ithread) * (scal(ist, ithread)/scal(ist, 1))**2
        end if
      end do
    end do

    ! the result is in scal(ist, 1) and ssq(ist, 1)
    do ist = 1, aa%nst
      nrm2(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        nrm2(ist) = hypot(nrm2(ist), scal(indb, 1)*sqrt(mesh%volume_element*ssq(indb, 1)))
      end do
    end do

    SAFE_DEALLOCATE_A(scal)
    SAFE_DEALLOCATE_A(ssq)

  case(BATCH_DEVICE_PACKED)

    ASSERT(.not. mesh%use_curvilinear)

    SAFE_ALLOCATE(ssq(1:aa%pack_size(1), 1))

    call accel_create_buffer(nrm2_buffer, ACCEL_MEM_WRITE_ONLY, TYPE_FLOAT, aa%pack_size(1))

    do ist = 1, aa%nst_linear
      call accel_set_stream(ist)
      call X(accel_nrm2)(N = int(mesh%np, 8), X = aa%ff_device, offx = int(ist - 1, 8), incx = int(aa%pack_size(1), 8), &
        res = nrm2_buffer, offres = int(ist - 1, 8))
    end do
    call accel_synchronize_all_streams()
    call accel_set_stream(1)

    call accel_read_buffer(nrm2_buffer, aa%pack_size(1), ssq)

    call accel_release_buffer(nrm2_buffer)

    do ist = 1, aa%nst
      nrm2(ist) = M_ZERO
      do idim = 1, aa%dim
        indb = aa%ist_idim_to_linear((/ist, idim/))
        nrm2(ist) = hypot(nrm2(ist), sqrt(mesh%volume_element)*ssq(indb, 1))
      end do
    end do

    SAFE_DEALLOCATE_A(ssq)

  end select
  
  ! REDUCTION IS REQUIRED, THIS IS DONE BY THE CALLING FUNCTION

  call profiling_out(prof)
  POP_SUB(X(priv_mesh_batch_nrm2))
end subroutine X(priv_mesh_batch_nrm2)

! ---------------------------------------------------------
!> Orthonormalizes states of phib to the orbitals of nst batches of psi.
!! It also permits doing only the orthogonalization (no normalization).
subroutine X(mesh_batch_orthogonalization)(mesh, nst, psib, phib,  &
  normalize, overlap, norm, gs_scheme)
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: nst
  class(batch_t),    intent(in)    :: psib(:)   !< psi(nst)
  class(batch_t),    intent(inout) :: phib      
  logical, optional, intent(in)    :: normalize
  R_TYPE,  optional, intent(out)   :: overlap(:,:) !< (nst, phib%nst)
  R_TYPE,  optional, intent(out)   :: norm(:)
  integer, optional, intent(in)    :: gs_scheme

  logical :: normalize_
  integer :: ist, is
  R_TYPE, allocatable   :: nrm2(:)
  R_TYPE, allocatable  :: ss(:,:), ss_full(:,:)
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof
  logical :: drcgs
  integer :: nsteps

  call profiling_in(prof, "BATCH_GRAM_SCHMIDT")
  PUSH_SUB(X(mesh_batch_orthogonalization))

  SAFE_ALLOCATE(ss(1:phib%nst, 1:nst))
  ss = R_TOTYPE(M_ZERO)

  do ist = 1, nst
    call phib%check_compatibility_with(psib(ist))
  end do

  drcgs = .false.
  nsteps = 1
  if(present(gs_scheme)) then
    if(gs_scheme == OPTION__ARNOLDIORTHOGONALIZATION__DRCGS) then
      drcgs = .true.
      nsteps = 2
      SAFE_ALLOCATE(ss_full(1:phib%nst, 1:nst))
      ss_full = R_TOTYPE(M_ZERO)
    end if
  end if

  do is = 1, nsteps
    if(nst>=1 .and. drcgs) then
      call X(mesh_batch_dotp_vector)(mesh, psib(nst), phib, ss(1:phib%nst,1))
      call batch_axpy(mesh%np, -ss(1:phib%nst,1), psib(nst), phib, a_full = .false.)
      if(present(overlap)) ss_full(1:phib%nst, nst) = ss_full(1:phib%nst, nst) + ss(1:phib%nst, 1)
    end if
    ss = R_TOTYPE(M_ZERO)

    !TODO: We should reuse phib here for improved performances
    do ist = 1, nst
      call X(mesh_batch_dotp_vector)(mesh, psib(ist), phib, ss(1:phib%nst,ist), reduce = .false.) 
    end do

    if(mesh%parallel_in_domains) then
      call profiling_in(reduce_prof, "BATCH_GRAM_SCHMIDT_REDUCE")
      call comm_allreduce(mesh%mpi_grp%comm, ss, dim = (/phib%nst, nst/))
      call profiling_out(reduce_prof)
    end if
   
    !TODO: We should have a routine batch_gemv for improved performances
    do ist = 1, nst
      call batch_axpy(mesh%np, -ss(1:phib%nst,ist), psib(ist), phib, a_full = .false.)
    end do

    !We accumulate the overlap
    if(drcgs .and. present(overlap)) then
      do ist = 1, nst
        ss_full(1:phib%nst, ist) = ss_full(1:phib%nst, ist) + ss(1:phib%nst, ist)
      end do 
    end if
  end do

  !We have a transpose here because this helps for the Lanczos implementation
  !which is the only routine using this one at the moment
  !Indeed, Lanczos acts on phib%nst arrays of dimension nst, whereas the code would return 
  !an array of dim (phib%nst, nst)
  !For an orthogalization, it is more natural to have for each state the overlap with the others
  !which is what the code outputs now.
  if(present(overlap)) then
    if(drcgs) then
      overlap(1:nst, 1:phib%nst) = transpose(ss_full(1:phib%nst, 1:nst))
    else
      overlap(1:nst, 1:phib%nst) = transpose(ss(1:phib%nst, 1:nst))
    end if
  end if

  normalize_ = optional_default(normalize, .false.)
  if(present(norm) .or. normalize_) then
    SAFE_ALLOCATE(nrm2(1:phib%nst))
    !Here we do not call mesh_batch_nrm2 which is too slow
    call X(mesh_batch_dotp_vector)(mesh, phib, phib, nrm2)
    if(present(norm)) then
      norm(1:phib%nst) = sqrt(TOFLOAT(nrm2(1:phib%nst)))
    end if
    if(normalize_) then
      call batch_scal(mesh%np, M_ONE/sqrt(TOFLOAT(nrm2)), phib, a_full =.false.)
    end if
    SAFE_DEALLOCATE_A(nrm2)
  end if

  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(ss_full)

  POP_SUB(X(mesh_batch_orthogonalization))
  call profiling_out(prof)
end subroutine X(mesh_batch_orthogonalization)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
