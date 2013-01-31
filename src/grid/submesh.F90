!! Copyright (C) 2007 X. Andrade
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
!! $Id: submesh.F90 2781 2007-03-23 10:58:32Z lorenzen $
 
#include "global.h"
  
module submesh_m
  use batch_m
  use blas_m
  use comm_m
  use global_m
  use lalg_basic_m
  use messages_m
  use mesh_m
  use mpi_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use solids_m
  use unit_m
  use unit_system_m
    
  implicit none
  private 

  public ::                      &
    submesh_t,                   &
    submesh_null,                &
    submesh_init_sphere,         &
    submesh_copy,                &
    submesh_get_inv,             &
    dsm_integrate,               &
    zsm_integrate,               &
    submesh_add_to_mesh,         &
    dsubmesh_batch_add,          &
    zsubmesh_batch_add,          &
    submesh_to_mesh_dotp,        &
    dsubmesh_batch_add_matrix,   &
    zsubmesh_batch_add_matrix,   &
    dsubmesh_batch_dotp_matrix,  &
    zsubmesh_batch_dotp_matrix,  &
    submesh_overlap,             &
    submesh_end

  type submesh_t
    FLOAT                 :: center(1:MAX_DIM)
    FLOAT                 :: radius
    integer               :: np = -1        !< number of points inside the submesh
    integer               :: np_part        !< number of points inside the submesh including ghost points
    integer,      pointer :: map(:) => null() !< index in the mesh of the points inside the sphere
    FLOAT,        pointer :: x(:,:) => null()
    type(mesh_t), pointer :: mesh
    logical               :: has_points
  end type submesh_t
  
  integer :: tagcounter = 0

  interface submesh_add_to_mesh
    module procedure ddsubmesh_add_to_mesh, zdsubmesh_add_to_mesh
  end interface submesh_add_to_mesh

  interface submesh_to_mesh_dotp
    module procedure ddsubmesh_to_mesh_dotp, zdsubmesh_to_mesh_dotp
  end interface submesh_to_mesh_dotp

contains
  
  subroutine submesh_null(sm)
    type(submesh_t), intent(out) :: sm

    PUSH_SUB(submesh_null)

    sm%np = -1
    nullify(sm%map)
    nullify(sm%x)
    nullify(sm%mesh)

    POP_SUB(submesh_null)

  end subroutine submesh_null

  ! -------------------------------------------------------------

  subroutine submesh_init_sphere(this, sb, mesh, center, rc)
    type(submesh_t),      intent(out)  :: this
    type(simul_box_t),    intent(in)   :: sb
    type(mesh_t), target, intent(in)   :: mesh
    FLOAT,                intent(in)   :: center(:)
    FLOAT,                intent(in)   :: rc
    
    FLOAT :: r2, xx(1:MAX_DIM)
    FLOAT, allocatable :: center_copies(:, :)
    integer :: icell, is, isb, ip, ix, iy, iz
    type(profile_t), save :: submesh_init_prof
    type(periodic_copy_t) :: pp
    integer, allocatable :: map_inv(:)
    integer :: nmax(1:MAX_DIM), nmin(1:MAX_DIM)

    PUSH_SUB(submesh_init_sphere)
    call profiling_in(submesh_init_prof, "SUBMESH_INIT")

    this%mesh => mesh

    this%center(1:sb%dim) = center(1:sb%dim)

    this%radius = rc

    ! The spheres are generated differently for periodic coordinates,
    ! mainly for performance reasons.
    if(.not. simul_box_is_periodic(sb)) then 

      SAFE_ALLOCATE(map_inv(0:this%mesh%np_part))
      map_inv(0:this%mesh%np_part) = 0
      
      nmin = 0
      nmax = 0

      ! get a cube of points that contains the sphere
      nmin(1:sb%dim) = int((center(1:sb%dim) - abs(rc) - sb%box_offset(1:sb%dim))/mesh%spacing(1:sb%dim)) - 1
      nmax(1:sb%dim) = int((center(1:sb%dim) + abs(rc) - sb%box_offset(1:sb%dim))/mesh%spacing(1:sb%dim)) + 1

      ! make sure that the cube is inside the grid
      nmin(1:sb%dim) = max(mesh%idx%nr(1, 1:sb%dim), nmin(1:sb%dim))
      nmax(1:sb%dim) = min(mesh%idx%nr(2, 1:sb%dim), nmax(1:sb%dim))

      ! Get the total number of points inside the sphere
      is = 0   ! this index counts inner points
      isb = 0  ! and this one boundary points
      do iz = nmin(3), nmax(3)
        do iy = nmin(2), nmax(2)
          do ix = nmin(1), nmax(1)
            ip = mesh%idx%lxyz_inv(ix, iy, iz)
#if defined(HAVE_MPI)
            if(ip == 0) cycle
            if(mesh%parallel_in_domains) ip = vec_global2local(mesh%vp, ip, mesh%vp%partno)
#endif
            if(ip == 0) cycle
            r2 = sum((mesh%x(ip, 1:sb%dim) - center(1:sb%dim))**2)
            if(r2 <= rc**2) then
              if(ip > mesh%np) then
                ! boundary points are marked as negative values
                isb = isb + 1
                map_inv(ip) = -isb
              else
                is = is + 1
                map_inv(ip) = is
              end if
            end if
          end do
        end do
      end do
      this%np = is
      this%np_part = is + isb
      
      SAFE_ALLOCATE(this%map(1:this%np_part))
      SAFE_ALLOCATE(this%x(1:this%np_part, 0:sb%dim))
      
      ! Generate the table and the positions
      do iz = nmin(3), nmax(3)
        do iy = nmin(2), nmax(2)
          do ix = nmin(1), nmax(1)
            ip = mesh%idx%lxyz_inv(ix, iy, iz)
#if defined(HAVE_MPI)
            if(ip == 0) cycle
            if(mesh%parallel_in_domains) ip = vec_global2local(mesh%vp, ip, mesh%vp%partno)
#endif
            is = map_inv(ip)
            if(is == 0) cycle
            if(is < 0) then
              ! it is a boundary point, move it to ns+1:ns_part range
              is = -is + this%np
              map_inv(ip) = is
            end if
            this%map(is) = ip
            this%x(is, 1:sb%dim) = mesh%x(ip, 1:sb%dim) - center(1:sb%dim)
            this%x(is, 0) = sqrt(sum(this%x(is, 1:sb%dim)**2))
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(map_inv)

    ! This is the case for a periodic system
    else

      ! Get the total number of points inside the sphere considering
      ! replicas along PBCs

      ! this requires some optimization

      if(any(rc >= sb%lsize(1:sb%periodic_dim))) then
        write(message(1), '(a,f5.2,3a)') 'Submesh radius (', &
          units_from_atomic(units_out%length, rc), " ", trim(units_abbrev(units_out%length)), &
          ") is larger than simulation box."
        write(message(2), '(a)') 'Multiple periodic replicas in the same region will not be considered.'
        call messages_warning(2)
      end if

      call periodic_copy_init(pp, sb, center(1:sb%dim), rc)
      
      SAFE_ALLOCATE(center_copies(1:sb%dim, 1:periodic_copy_num(pp)))

      do icell = 1, periodic_copy_num(pp)
        center_copies(1:sb%dim, icell) = periodic_copy_position(pp, sb, icell)
      end do

      is = 0
      do ip = 1, mesh%np_part
        do icell = 1, periodic_copy_num(pp)
          r2 = sum((mesh%x(ip, 1:sb%dim) - center_copies(1:sb%dim, icell))**2)
          if(r2 > rc**2 ) cycle
          is = is + 1
        end do
        if (ip == mesh%np) this%np = is
      end do
      
      this%np_part = is

      SAFE_ALLOCATE(this%map(1:this%np_part))
      SAFE_ALLOCATE(this%x(1:this%np_part, 0:sb%dim))
            
      !iterate again to fill the tables
      is = 0
      do ip = 1, mesh%np_part
        do icell = 1, periodic_copy_num(pp)
          xx(1:sb%dim) = mesh%x(ip, 1:sb%dim) - center_copies(1:sb%dim, icell)
          r2 = sum(xx(1:sb%dim)**2)
          if(r2 > rc**2 ) cycle
          is = is + 1
          this%map(is) = ip
          this%x(is, 0) = sqrt(r2)
          this%x(is, 1:sb%dim) = xx(1:sb%dim)
         end do
      end do

      SAFE_DEALLOCATE_A(center_copies)
      
      call periodic_copy_end(pp)

    end if

    this%has_points = (this%np > 0)

    call profiling_out(submesh_init_prof)
    POP_SUB(submesh_init_sphere)
  end subroutine submesh_init_sphere

  ! --------------------------------------------------------------

  subroutine submesh_end(this)
    type(submesh_t),   intent(inout)  :: this
    
    PUSH_SUB(submesh_end)

    if( this%np /= -1 ) then
      nullify(this%mesh)
      this%np = -1
      SAFE_DEALLOCATE_P(this%map)
      SAFE_DEALLOCATE_P(this%x)
    end if

    POP_SUB(submesh_end)

  end subroutine submesh_end

  ! --------------------------------------------------------------

  subroutine submesh_copy(sm_in, sm_out)
    type(submesh_t), target,  intent(in)   :: sm_in
    type(submesh_t),          intent(out)  :: sm_out

    PUSH_SUB(submesh_copy)
    
    ASSERT(sm_out%np == -1)

    sm_out%mesh => sm_in%mesh

    sm_out%center = sm_in%center
    sm_out%radius = sm_in%radius

    sm_out%np = sm_in%np
    sm_out%np_part = sm_in%np_part
    
    SAFE_ALLOCATE(sm_out%map(1:sm_out%np_part))
    SAFE_ALLOCATE(sm_out%x(1:sm_out%np_part, 0:ubound(sm_in%x, 2)))

    sm_out%map(1:sm_out%np_part) = sm_in%map(1:sm_in%np_part)
    sm_out%x(1:sm_out%np_part, 0:ubound(sm_in%x, 2)) = sm_in%x(1:sm_in%np_part, 0:ubound(sm_in%x, 2))

    POP_SUB(submesh_copy)

  end subroutine submesh_copy

  ! --------------------------------------------------------------

  subroutine submesh_get_inv(this, map_inv)
    type(submesh_t),      intent(in)   :: this
    integer,              intent(out)  :: map_inv(:)

    integer :: is

    PUSH_SUB(submesh_get_inv)
    
    map_inv(1:this%mesh%np_part) = 0
    forall (is = 1:this%np) map_inv(this%map(is)) = is

    POP_SUB(submesh_get_inv)
  end subroutine submesh_get_inv

  ! --------------------------------------------------------------

  logical pure function submesh_overlap(sm1, sm2) result(overlap)
    type(submesh_t),      intent(in)   :: sm1
    type(submesh_t),      intent(in)   :: sm2
    
    FLOAT :: distance

    distance = sum((sm1%center(1:MAX_DIM) - sm2%center(1:MAX_DIM))**2)
    overlap = distance + CNST(100.0)*M_EPSILON <= (sm1%radius + sm2%radius)**2

  end function submesh_overlap
  


#include "undef.F90"
#include "real.F90"
#include "submesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "submesh_inc.F90"

end module submesh_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
