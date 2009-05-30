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
    SAFE_ALLOCATE(ddtmp(1:aa%nst, 1:bb%nst))
    forall(ist = 1:aa%nst, jst = 1:bb%nst) ddtmp(ist, jst) = dd(ist, jst)
    call MPI_Allreduce(ddtmp, dd, aa%nst*bb%nst, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(ddtmp)
    call profiling_out(profcomm)
  end if
#endif

  forall(ist = 1:aa%nst, jst = 1:bb%nst) dot(aa%states(ist)%ist, bb%states(jst)%ist) = dd(ist, jst)

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  call pop_sub()
end subroutine X(mesh_batch_dotp_matrix)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
