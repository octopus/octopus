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

  integer :: ist, jst, idim, sp, block_size, ep, ip, lda, ldb
  R_TYPE :: ss
  R_TYPE, allocatable :: dd(:, :)
#ifdef HAVE_MPI
  R_TYPE, allocatable :: ddtmp(:, :)
#endif
  type(profile_t), save :: prof, profgemm, profcomm
  logical :: use_blas, reduce_

  call push_sub('mesh_function_inc.Xmf_dotp_batch')
  call profiling_in(prof, "DOTP_BATCH")

#ifdef HAVE_MPI
  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce
#endif

  ASSERT(aa%dim == bb%dim)

  SAFE_ALLOCATE(dd(1:aa%nst, 1:bb%nst))

  use_blas = associated(aa%X(psicont)) .and. associated(bb%X(psicont)) .and. (.not. mesh%use_curvilinear) .and. (aa%dim == 1)

  if(use_blas) then
    call profiling_in(profgemm, "DOTP_BATCH_GEMM")

    lda = size(aa%X(psicont), dim = 1)
    ldb = size(bb%X(psicont), dim = 1)
    call blas_gemm('c', 'n', aa%nst, bb%nst, mesh%np, &
         R_TOTYPE(mesh%vol_pp(1)), &
         aa%X(psicont)(1, 1, 1), lda, &
         bb%X(psicont)(1, 1, 1), ldb, &
         R_TOTYPE(M_ZERO), dd(1, 1), aa%nst)

  else

    dd = R_TOTYPE(M_ZERO)

    block_size = hardware%X(block_size)

    do idim = 1, aa%dim
      do sp = 1, mesh%np, block_size
        ep = min(mesh%np, sp + block_size - 1)

        if(mesh%use_curvilinear) then

          do ist = 1, aa%nst
            do jst = 1, bb%nst

              ss = M_ZERO
              do ip = sp, ep
                ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%states(ist)%X(psi)(ip, idim))*bb%states(jst)%X(psi)(ip, idim)
              end do
              dd(ist, jst) = dd(ist, jst) + ss

            end do
          end do

        else

          do ist = 1, aa%nst
            do jst = 1, bb%nst
              dd(ist, jst) = dd(ist, jst) + &
                   mesh%vol_pp(1)*blas_dot(ep - sp + 1, aa%states(ist)%X(psi)(sp, idim), 1, bb%states(jst)%X(psi)(sp, idim), 1)
            end do
          end do

        end if
      end do
    end do

  end if

  if(mesh%use_curvilinear) then
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + 2*R_MUL))
  else
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + R_MUL))
  end if

  if(use_blas) call profiling_out(profgemm)

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "DOTP_BATCH_REDUCE")
#ifndef HAVE_MPI2
    SAFE_ALLOCATE(ddtmp(1:aa%nst, 1:bb%nst))
    forall(ist = 1:aa%nst, jst = 1:bb%nst) ddtmp(ist, jst) = dd(ist, jst)
#endif
    call MPI_Allreduce(MPI_IN_PLACE_OR(ddtmp), dd, aa%nst*bb%nst, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(ddtmp)
    call profiling_out(profcomm)
  end if
#endif

  forall(ist = 1:aa%nst, jst = 1:bb%nst) dot(aa%states(ist)%ist, bb%states(jst)%ist) = dd(ist, jst)

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  call pop_sub('mesh_function_inc.Xmf_dotp_batch')
end subroutine X(mesh_batch_dotp_matrix)

!-----------------------------------------------------------------

subroutine X(mesh_batch_dotp_self)(mesh, aa, dot, reduce)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(in)    :: aa
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: reduce

  integer :: ist, jst, idim, sp, block_size, ep, ip, lda
  R_TYPE :: ss
#ifdef HAVE_MPI
  R_TYPE, allocatable :: dottmp(:, :)
#endif
  type(profile_t), save :: prof, profgemm, profcomm
  logical :: use_blas, reduce_

  call push_sub('mesh_function_inc.Xmf_dotp_batch')
  call profiling_in(prof, "BATCH_DOTP_SELF")

  ! some limitations of the current implementation
  ASSERT(ubound(dot, dim = 1) == aa%nst .and. ubound(dot, dim = 2) == aa%nst)
  ASSERT(aa%states(1)%ist == 1)

#ifdef HAVE_MPI
  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce
#endif

  use_blas = associated(aa%X(psicont)) .and. (.not. mesh%use_curvilinear)

  if(use_blas) then
    call profiling_in(profgemm, "BATCH_HERK")

    lda = size(aa%X(psicont), dim = 1)*aa%dim

    call blas_herk('l', 'c', aa%nst, mesh%np, R_TOTYPE(mesh%vol_pp(1)), aa%X(psicont)(1, 1, 1), &
      lda, R_TOTYPE(M_ZERO), dot(1, 1), aa%nst)

    if(aa%dim == 2) then
      call blas_herk('l', 'c', aa%nst, mesh%np, R_TOTYPE(mesh%vol_pp(1)), aa%X(psicont)(1, 2, 1), &
        lda, R_TOTYPE(M_ONE), dot(1, 1), aa%nst)
    end if

  else

    dot = R_TOTYPE(M_ZERO)

    block_size = hardware%X(block_size)

    do idim = 1, aa%dim
      do sp = 1, mesh%np, block_size
        ep = min(mesh%np, sp + block_size - 1)

        if(mesh%use_curvilinear) then

          do ist = 1, aa%nst
            do jst = 1, ist

              ss = M_ZERO
              do ip = sp, ep
                ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%states(ist)%X(psi)(ip, idim))*aa%states(jst)%X(psi)(ip, idim)
              end do
              dot(ist, jst) = dot(ist, jst) + ss

            end do
          end do

        else

          do ist = 1, aa%nst
            do jst = 1, ist
              dot(ist, jst) = dot(ist, jst) + &
                mesh%vol_pp(1)*blas_dot(ep - sp + 1, aa%states(ist)%X(psi)(sp, idim), 1, aa%states(jst)%X(psi)(sp, idim), 1)
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

  forall(ist = 1:aa%nst)
    forall(jst = 1:ist - 1) dot(jst, ist) = R_CONJ(dot(ist, jst))
  end forall

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains .and. reduce_) then
    call profiling_in(profcomm, "BATCH_SELF_REDUCE")
#ifndef HAVE_MPI2
    SAFE_ALLOCATE(dottmp(1:aa%nst, 1:aa%nst))
    forall(ist = 1:aa%nst, jst = 1:aa%nst) dottmp(ist, jst) = dot(ist, jst)
#endif
    call MPI_Allreduce(MPI_IN_PLACE_OR(dottmp), dot, aa%nst**2, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(dottmp)
    call profiling_out(profcomm)
  end if
#endif

  call profiling_out(prof)
  call pop_sub('mesh_function_inc.Xmf_dotp_batch')
end subroutine X(mesh_batch_dotp_self)

!-----------------------------------------------------------------

subroutine X(mesh_batch_rotate)(mesh, aa, transf)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(inout) :: aa
  R_TYPE,            intent(inout) :: transf(:, :)

  R_TYPE, allocatable :: psinew(:, :), psicopy(:, :)
  
  integer :: ist, idim, block_size, size, sp
  type(profile_t), save :: prof
  logical :: use_blas

  call profiling_in(prof, "BATCH_ROTATE")

#ifdef R_TREAL  
  block_size = max(40, hardware%l2%size/(2*8*aa%nst))
#else
  block_size = max(20, hardware%l2%size/(2*16*aa%nst))
#endif

  use_blas = associated(aa%X(psicont))

  SAFE_ALLOCATE(psinew(1:block_size, 1:aa%nst))

  if(.not. use_blas) then
    SAFE_ALLOCATE(psicopy(1:block_size, 1:aa%nst))
  end if

  do sp = 1, mesh%np, block_size
    size = min(block_size, mesh%np - sp + 1)
    
    do idim = 1, aa%dim
      
      if(use_blas) then
        
        call blas_gemm('N', 'N', &
          size, aa%nst, aa%nst, &
          R_TOTYPE(M_ONE), aa%X(psicont)(sp, idim, 1), ubound(aa%X(psicont), dim=1)*aa%dim, &
          transf(1, 1), aa%nst, &
          R_TOTYPE(M_ZERO), psinew(1, 1), block_size)
        
        do ist = 1, aa%nst
          call blas_copy(size, psinew(1, ist), 1, aa%X(psicont)(sp, idim, ist), 1)
        end do
        
      else

        do ist = 1, aa%nst
          call blas_copy(size, aa%states(ist)%X(psi)(sp, idim), 1, psicopy(1, ist), 1)
        end do

        call blas_gemm('N', 'N', &
          size, aa%nst, aa%nst, &
          R_TOTYPE(M_ONE), psicopy(1, 1), block_size, &
          transf(1, 1), aa%nst, &
          R_TOTYPE(M_ZERO), psinew(1, 1), block_size)

        do ist = 1, aa%nst
          call blas_copy(size, psinew(1, ist), 1, aa%states(ist)%X(psi)(sp, idim), 1)
        end do

      end if
    end do
  end do

  SAFE_DEALLOCATE_A(psicopy)
  SAFE_DEALLOCATE_A(psinew)

  call profiling_count_operations((R_ADD + R_MUL)*dble(mesh%np)*aa%dim*aa%nst**2)

  call profiling_out(prof)

end subroutine X(mesh_batch_rotate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
