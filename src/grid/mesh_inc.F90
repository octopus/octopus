!! Copyright (C) 2020 S. Ohlmann
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

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_0)(mesh, aa)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa

  call comm_allreduce(mesh%mpi_grp, aa)
end subroutine X(mesh_allreduce_0)

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_1)(mesh, aa, dim)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa(:)
  integer, optional,              intent(in)    :: dim

  call comm_allreduce(mesh%mpi_grp, aa, dim)
end subroutine X(mesh_allreduce_1)

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_2)(mesh, aa, dim)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa(:, :)
  integer, optional,              intent(in)    :: dim(:) !< (2)

  call comm_allreduce(mesh%mpi_grp, aa, dim)
end subroutine X(mesh_allreduce_2)

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_3)(mesh, aa, dim)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa(:, :, :)
  integer, optional,              intent(in)    :: dim(:) !< (3)

  call comm_allreduce(mesh%mpi_grp, aa, dim)
end subroutine X(mesh_allreduce_3)

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_4)(mesh, aa)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa(:, :, :, :)

  call comm_allreduce(mesh%mpi_grp, aa)
end subroutine X(mesh_allreduce_4)

! -----------------------------------------------------------------------------
subroutine X(mesh_allreduce_5)(mesh, aa)
  class(mesh_t),                  intent(in)    :: mesh
  R_TYPE,                         intent(inout) :: aa(:, :, :, :, :)

  call comm_allreduce(mesh%mpi_grp, aa)
end subroutine X(mesh_allreduce_5)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
