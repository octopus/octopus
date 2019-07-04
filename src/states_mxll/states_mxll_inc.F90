!! Copyright (C) 2011 X. Andrade
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
!----------------------------------------------------------


!----------------------------------------------------------
subroutine X(ma_mesh_to_mx_mesh)(maxwell_st, maxwell_gr, st, gr, ma_field_in, mx_field_out, field_dim)
  type(states_mxll_t),    intent(in)    :: maxwell_st
  type(grid_t),           intent(in)    :: maxwell_gr
  type(states_t),         intent(in)    :: st
  type(grid_t),           intent(in)    :: gr
  R_TYPE,                 intent(inout) :: ma_field_in(:,:)
  R_TYPE,                 intent(inout) :: mx_field_out(:,:)
  integer,                intent(in)    :: field_dim

  integer             :: idim, idx, ip_ma_local, ip_ma_global, ip_mx_local, ip_mx_global, mpi_err
  R_TYPE, allocatable :: ma_field_in_global(:,:)

  PUSH_SUB(X(ma_mesh_to_mx_mesh))

  SAFE_ALLOCATE(ma_field_in_global(1:gr%mesh%np_part_global,field_dim))

  ma_field_in_global(:,:) = M_z0
  if (gr%mesh%parallel_in_domains) then
    do idim=1, field_dim
#if defined(HAVE_MPI)
      call vec_allgather(gr%mesh%vp, ma_field_in_global(:,idim), ma_field_in(:,idim))
      call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
  else
    ma_field_in_global(:,:) = ma_field_in(:,:)
  end if

  mx_field_out(:,:) = M_z0
  if (gr%mesh%parallel_in_domains) then
    if (maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap_number > 0) then
      do idx=1, maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap_number
        ip_mx_local  = maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap(idx)
        ip_mx_global = vec_local2global_part(maxwell_gr%mesh%vp, ip_mx_local, maxwell_gr%mesh%vp%partno)
        ip_ma_global = maxwell_gr%mesh%mx_ma_mesh_mapping%global2_overlap(idx)
        if ((gr%mesh%idx%lxyz(ip_ma_global,1) /= maxwell_gr%mesh%idx%lxyz(ip_mx_global,1)) .or. &
            (gr%mesh%idx%lxyz(ip_ma_global,2) /= maxwell_gr%mesh%idx%lxyz(ip_mx_global,2)) .or. &
            (gr%mesh%idx%lxyz(ip_ma_global,3) /= maxwell_gr%mesh%idx%lxyz(ip_mx_global,3))) then
          message(1) = "There is an error in the Maxwell grid point to the matter grid point!"
        end if
        mx_field_out(ip_mx_local,:) = ma_field_in_global(ip_ma_global,:)
      end do
    end if
  else
    if (maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap_number > 0) then
      do idx=1, maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap_number
        ip_mx_local = maxwell_gr%mesh%mx_ma_mesh_mapping%local1_overlap(idx)
        ip_ma_local = maxwell_gr%mesh%mx_ma_mesh_mapping%global2_overlap(idx)
        mx_field_out(ip_mx_local,:) = ma_field_in(ip_ma_local,:)
      end do
    end if
  end if

  SAFE_DEALLOCATE_A(ma_field_in_global)

  POP_SUB(X(ma_mesh_to_mx_mesh))
end subroutine X(ma_mesh_to_mx_mesh)



!----------------------------------------------------------
subroutine X(mx_mesh_to_ma_mesh)(maxwell_st, maxwell_gr, st, gr, mx_field_in, ma_field_out, field_dim)
  type(states_mxll_t),      intent(in)    :: maxwell_st
  type(grid_t),        intent(in)    :: maxwell_gr
  type(states_t),      intent(in)    :: st
  type(grid_t),        intent(in)    :: gr
  R_TYPE,               intent(inout) :: mx_field_in(:,:)
  R_TYPE,               intent(inout) :: ma_field_out(:,:)
  integer,             intent(in)    :: field_dim

  integer            :: idim, idx, ip_mx_local, ip_mx_global, ip_ma_local, ip_ma_global, mpi_err
  R_TYPE, allocatable :: mx_field_in_global(:,:)

  PUSH_SUB(X(mx_mesh_to_ma_mesh))

  SAFE_ALLOCATE(mx_field_in_global(1:maxwell_gr%mesh%np_part_global,field_dim))

  mx_field_in_global(:,:) = M_z0
  if (maxwell_gr%mesh%parallel_in_domains) then
    do idim=1, field_dim
#if defined(HAVE_MPI)
      call vec_allgather(maxwell_gr%mesh%vp, mx_field_in_global(:,idim), mx_field_in(:,idim))
      call MPI_Barrier(maxwell_gr%mesh%vp%comm, mpi_err)
#endif
      end do
  else
    mx_field_in_global(:,:) = mx_field_in(:,:)
  end if

  ma_field_out(:,:) = M_z0
  if (maxwell_gr%mesh%parallel_in_domains) then
    if (gr%mesh%ma_mx_mesh_mapping%local1_overlap_number > 0) then
      do idx=1, gr%mesh%ma_mx_mesh_mapping%local1_overlap_number
        ip_ma_local  = gr%mesh%ma_mx_mesh_mapping%local1_overlap(idx)
        ip_ma_global = vec_local2global_part(gr%mesh%vp, ip_ma_local, maxwell_gr%mesh%vp%partno)
        ip_mx_global = gr%mesh%ma_mx_mesh_mapping%global2_overlap(idx)
        if ((maxwell_gr%mesh%idx%lxyz(ip_mx_global,1) /= gr%mesh%idx%lxyz(ip_ma_global,1)) .or. &
            (maxwell_gr%mesh%idx%lxyz(ip_mx_global,2) /= gr%mesh%idx%lxyz(ip_ma_global,2)) .or. &
            (maxwell_gr%mesh%idx%lxyz(ip_mx_global,3) /= gr%mesh%idx%lxyz(ip_ma_global,3))) then
          message(1) = "There is an error in the Maxwell grid point to the matter grid point!"
        end if
        ma_field_out(ip_ma_local,:) = mx_field_in_global(ip_mx_global,:)
      end do
    end if
  else
    if (gr%mesh%ma_mx_mesh_mapping%local1_overlap_number > 0) then
      do idx=1, gr%mesh%ma_mx_mesh_mapping%local1_overlap_number
        ip_ma_local = gr%mesh%ma_mx_mesh_mapping%local1_overlap(idx)
        ip_mx_local = gr%mesh%ma_mx_mesh_mapping%global2_overlap(idx)
        ma_field_out(ip_ma_local,:) = mx_field_in(ip_mx_local,:)
      end do
    end if
  end if

  SAFE_DEALLOCATE_A(mx_field_in_global)

  POP_SUB(X(mx_mesh_to_ma_mesh))
end subroutine X(mx_mesh_to_ma_mesh)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
