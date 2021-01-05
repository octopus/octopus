!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module symmetrizer_oct_m
  use comm_oct_m
  use global_oct_m
  use index_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use symm_op_oct_m
  use symmetries_oct_m

  implicit none

  private
  public ::                             &
    symmetrizer_t,                      &
    symmetrizer_init,                   &
    symmetrizer_end,                    &
    dsymmetrizer_apply,                 &
    zsymmetrizer_apply,                 &
    dsymmetrizer_apply_single,          &
    zsymmetrizer_apply_single,          &
    dsymmetrize_tensor_cart,            &
    zsymmetrize_tensor_cart,            &
    dsymmetrize_magneto_optics_cart,    &
    zsymmetrize_magneto_optics_cart

  type symmetrizer_t
    private
    type(mesh_t), pointer :: mesh
    integer, allocatable :: map(:,:)
    integer, allocatable :: map_inv(:,:)
  end type symmetrizer_t

contains

  ! ---------------------------------------------------------
  subroutine symmetrizer_init(this, mesh)
    type(symmetrizer_t),         intent(out) :: this
    type(mesh_t),        target, intent(in)  :: mesh

    integer :: nops, ip, iop, idir, idx(MAX_DIM)
    FLOAT :: destpoint(1:3), srcpoint(1:3), srcpoint_inv(1:3), lsize(1:3), offset(1:3)
    type(profile_t), save :: prof

    PUSH_SUB(symmetrizer_init)
    
    this%mesh => mesh

    !For each operation, we create a mapping between the grid point and the symmetric point
    nops = symmetries_number(mesh%sb%symm)

    SAFE_ALLOCATE(this%map(1:mesh%np, 1:nops))
    SAFE_ALLOCATE(this%map_inv(1:mesh%np, 1:nops))

    call profiling_in(prof, "SYMMETRIZER_INIT") 

    lsize(1:3) = TOFLOAT(mesh%idx%ll(1:3))
    offset(1:3) = TOFLOAT(mesh%idx%nr(1, 1:3) + mesh%idx%enlarge(1:3))

    do ip = 1, mesh%np
      if(mesh%parallel_in_domains) then
        ! convert to global point
        call index_to_coords(mesh%idx, mesh%vp%local(mesh%vp%xlocal + ip - 1), idx)
      else
        call index_to_coords(mesh%idx, ip, idx)
      end if
      destpoint(1:3) = TOFLOAT(idx(1:3)) - offset(1:3)
      ! offset moves corner of cell to origin, in integer mesh coordinates

      ASSERT(all(destpoint >= 0))
      ASSERT(all(destpoint < lsize))

      ! move to center of cell in real coordinates
      destpoint = destpoint - TOFLOAT(int(lsize)/2)

      !convert to proper reduced coordinates
      do idir = 1, 3
        destpoint(idir) = destpoint(idir)/lsize(idir)
      end do

      ! iterate over all points that go to this point by a symmetry operation
      do iop = 1, nops
        srcpoint = symm_op_apply_red(mesh%sb%symm%ops(iop), destpoint)
        srcpoint_inv = symm_op_apply_inv_red(mesh%sb%symm%ops(iop), destpoint)

        !We now come back to what should be an integer, if the symmetric point beloings to the grid
        !At this point, this is already checked
        do idir = 1, 3
          srcpoint(idir) = srcpoint(idir)*lsize(idir)
          srcpoint_inv(idir) = srcpoint_inv(idir)*lsize(idir)
        end do

        ! move back to reference to origin at corner of cell
        srcpoint = srcpoint + TOFLOAT(int(lsize)/2)
        srcpoint_inv = srcpoint_inv + TOFLOAT(int(lsize)/2)

        ! apply periodic boundary conditions in periodic directions
        do idir = 1, mesh%sb%periodic_dim
          if(srcpoint(idir) < M_ZERO .or. srcpoint(idir) + M_HALF*SYMPREC >= lsize(idir)) then
            srcpoint(idir) = modulo(srcpoint(idir)+M_HALF*SYMPREC, lsize(idir))
          end if
          if(srcpoint_inv(idir) < M_ZERO .or. srcpoint_inv(idir) + M_HALF*SYMPREC >= lsize(idir)) then
            srcpoint_inv(idir) = modulo(srcpoint_inv(idir)+M_HALF*SYMPREC, lsize(idir))
          end if
        end do
        ASSERT(all(srcpoint >= -SYMPREC))
        ASSERT(all(srcpoint < lsize))
        srcpoint(1:3) = srcpoint(1:3) + offset(1:3)

        ASSERT(all(srcpoint_inv >= -SYMPREC))
        ASSERT(all(srcpoint_inv < lsize))
        srcpoint_inv(1:3) = srcpoint_inv(1:3) + offset(1:3)

        this%map(ip, iop) = index_from_coords(this%mesh%idx, &
          [nint(srcpoint(1)), nint(srcpoint(2)), nint(srcpoint(3))])
        this%map_inv(ip, iop) = index_from_coords(this%mesh%idx, &
          [nint(srcpoint_inv(1)), nint(srcpoint_inv(2)), nint(srcpoint_inv(3))])
      end do
    end do

    call profiling_out(prof)

    POP_SUB(symmetrizer_init)
  end subroutine symmetrizer_init

  ! ---------------------------------------------------------

  subroutine symmetrizer_end(this)
    type(symmetrizer_t), intent(inout) :: this

    PUSH_SUB(symmetrizer_end)
    nullify(this%mesh)

    SAFE_DEALLOCATE_A(this%map)
    SAFE_DEALLOCATE_A(this%map_inv)

    POP_SUB(symmetrizer_end)
  end subroutine symmetrizer_end

  ! ---------------------------------------------------------
  
#include "undef.F90"
#include "real.F90"
#include "symmetrizer_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "symmetrizer_inc.F90"

end module symmetrizer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
