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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
 
#include "global.h"
  
module submesh_oct_m
  use accel_oct_m
  use batch_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use messages_oct_m
  use sort_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use types_oct_m
    
  implicit none
  private 

  public ::                      &
    submesh_t,                   &
    submesh_null,                &
    submesh_init,                &
    submesh_merge,               &
    submesh_shift_center,        &
    submesh_broadcast,           &    
    submesh_copy,                &
    submesh_get_inv,             &
    submesh_build_global,        &
    submesh_end_global,          &
    dsm_integrate,               &
    zsm_integrate,               &
    dsm_integrate_frommesh,      &
    zsm_integrate_frommesh,      & 
    dsm_nrm2,                    &
    zsm_nrm2,                    &
    submesh_add_to_mesh,         &
    dsubmesh_batch_add,          &
    zsubmesh_batch_add,          &
    submesh_to_mesh_dotp,        &
    dsubmesh_batch_add_matrix,   &
    zsubmesh_batch_add_matrix,   &
    dsubmesh_batch_dotp_matrix,  &
    zsubmesh_batch_dotp_matrix,  &
    submesh_overlap,             &
    dsubmesh_copy_from_mesh,     &
    zsubmesh_copy_from_mesh,     &
    dsubmesh_copy_from_mesh_batch,     &
    zsubmesh_copy_from_mesh_batch,     &
    submesh_end,                 &
    submesh_get_cube_dim,        &
    submesh_init_cube_map,       &
    submesh_end_cube_map

  type submesh_t
    ! Components are public by default
    FLOAT                 :: center(1:MAX_DIM)
    FLOAT                 :: radius
    integer               :: np             !< number of points inside the submesh
    integer               :: np_part        !< number of points inside the submesh including ghost points
    integer,      pointer :: map(:)         !< index in the mesh of the points inside the sphere
    FLOAT,        pointer :: x(:,:)
    type(mesh_t), pointer :: mesh
    logical               :: overlap        !< .true. if the submesh has more than one point that is mapped to a mesh point
    integer               :: np_global      !< total number of points in the entire mesh
    FLOAT,    allocatable :: x_global(:,:)  
    integer,  allocatable :: part_v(:)
    integer,  allocatable :: global2local(:)
    
    type(mesh_cube_map_t) :: cube_map
  end type submesh_t
  
  interface submesh_add_to_mesh
    module procedure ddsubmesh_add_to_mesh, zdsubmesh_add_to_mesh, zzsubmesh_add_to_mesh
  end interface submesh_add_to_mesh

  interface submesh_to_mesh_dotp
    module procedure ddsubmesh_to_mesh_dotp, zdsubmesh_to_mesh_dotp, zzsubmesh_to_mesh_dotp
  end interface submesh_to_mesh_dotp

   type(profile_t), save ::           &
       C_PROFILING_SM_REDUCE,         &
       C_PROFILING_SM_NRM2

contains
  
  subroutine submesh_null(sm)
    type(submesh_t), intent(inout) :: sm !< valgrind objects to intent(out) due to the initializations above

    PUSH_SUB(submesh_null)

    sm%np = -1
    sm%radius = M_ZERO
    nullify(sm%map)
    nullify(sm%x)
    nullify(sm%mesh)

    sm%np_global = -1

    POP_SUB(submesh_null)

  end subroutine submesh_null

  ! -------------------------------------------------------------

  subroutine submesh_init(this, sb, mesh, center, rc)
    type(submesh_t),      intent(inout)  :: this !< valgrind objects to intent(out) due to the initializations above
    type(simul_box_t),    intent(in)     :: sb
    type(mesh_t), target, intent(in)     :: mesh
    FLOAT,                intent(in)     :: center(:)
    FLOAT,                intent(in)     :: rc
    
    FLOAT :: r2, xx(1:MAX_DIM)
    FLOAT, allocatable :: center_copies(:, :), xtmp(:, :)
    integer :: icell, is, isb, ip, ix, iy, iz
    type(profile_t), save :: submesh_init_prof
    type(periodic_copy_t) :: pp
    integer, allocatable :: map_inv(:)
    integer :: nmax(1:MAX_DIM), nmin(1:MAX_DIM)
    integer, allocatable :: order(:)

    
    PUSH_SUB(submesh_init)
    call profiling_in(submesh_init_prof, "SUBMESH_INIT")

    call submesh_null(this)

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
      nmin(1:sb%dim) = int((center(1:sb%dim) - abs(rc))/mesh%spacing(1:sb%dim)) - 1
      nmax(1:sb%dim) = int((center(1:sb%dim) + abs(rc))/mesh%spacing(1:sb%dim)) + 1

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
      SAFE_ALLOCATE(xtmp(1:this%np_part, 0:sb%dim))
      
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
            xtmp(is, 1:sb%dim) = mesh%x(ip, 1:sb%dim) - center(1:sb%dim)
            xtmp(is, 0) = sqrt(sum(xtmp(is, 1:sb%dim)**2))
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(map_inv)

    ! This is the case for a periodic system
    else

      ! Get the total number of points inside the sphere considering
      ! replicas along PBCs

      ! this requires some optimization

      call periodic_copy_init(pp, sb, center(1:sb%dim), rc)
      
      SAFE_ALLOCATE(center_copies(1:sb%dim, 1:periodic_copy_num(pp)))

      do icell = 1, periodic_copy_num(pp)
        center_copies(1:sb%dim, icell) = periodic_copy_position(pp, sb, icell)
      end do

      is = 0
      do ip = 1, mesh%np_part
        do icell = 1, periodic_copy_num(pp)
          xx(1:sb%dim) = mesh%x(ip, 1:sb%dim) - center_copies(1:sb%dim, icell)
          r2 = sum(xx(1:sb%dim)**2)
          if(r2 > rc**2 ) cycle
          is = is + 1
        end do
        if (ip == mesh%np) this%np = is
      end do
      
      this%np_part = is

      SAFE_ALLOCATE(this%map(1:this%np_part))
      SAFE_ALLOCATE(xtmp(1:this%np_part, 0:sb%dim))
            
      !iterate again to fill the tables
      is = 0
      do ip = 1, mesh%np_part
        do icell = 1, periodic_copy_num(pp)
          xx(1:sb%dim) = mesh%x(ip, 1:sb%dim) - center_copies(1:sb%dim, icell)
          r2 = sum(xx(1:sb%dim)**2)
          if(r2 > rc**2 ) cycle
          is = is + 1
          this%map(is) = ip
          xtmp(is, 0) = sqrt(r2)
          xtmp(is, 1:sb%dim) = xx(1:sb%dim)
         end do
      end do

      SAFE_DEALLOCATE_A(center_copies)
      
      call periodic_copy_end(pp)

    end if

    ! now order points for better locality
    
    SAFE_ALLOCATE(order(1:this%np_part))
    SAFE_ALLOCATE(this%x(1:this%np_part, 0:sb%dim))

    forall(ip = 1:this%np_part) order(ip) = ip

    call sort(this%map, order)

    forall(ip = 1:this%np_part) this%x(ip, 0:sb%dim) = xtmp(order(ip), 0:sb%dim)

    !check whether points overlap
    
    this%overlap = .false.
    do ip = 1, this%np_part - 1
      if(this%map(ip) == this%map(ip + 1)) then
        this%overlap = .true.
        exit
      end if
    end do

    SAFE_DEALLOCATE_A(order)
    SAFE_DEALLOCATE_A(xtmp)

    call profiling_out(submesh_init_prof)
    POP_SUB(submesh_init)
  end subroutine submesh_init

  ! --------------------------------------------------------------
  !This routine takes two submeshes and merge them into a bigger submesh
  !The grid is centered on the first center
  subroutine submesh_merge(this, sb, mesh, sm1, sm2, shift)
    type(submesh_t),      intent(inout)  :: this !< valgrind objects to intent(out) due to the initializations above
    type(simul_box_t),    intent(in)     :: sb
    type(mesh_t), target, intent(in)     :: mesh
    type(submesh_t),      intent(in)     :: sm1
    type(submesh_t),      intent(in)     :: sm2
    FLOAT, optional,      intent(in)     :: shift(:) !< If present, shifts the center of sm2
    
    FLOAT :: r2 
    integer :: ip, is
    type(profile_t), save :: prof
    FLOAT :: xx(1:MAX_DIM), diff_centers(1:MAX_DIM)
    
    PUSH_SUB(submesh_merge)
    call profiling_in(prof, "SUBMESH_MERGE")

    this%mesh => mesh

    this%center(1:sb%dim)  = sm1%center(1:sb%dim)
    this%radius = sm1%radius

    diff_centers(1:sb%dim) = sm1%center(1:sb%dim)-sm2%center(1:sb%dim)
    if(present(shift)) diff_centers(1:sb%dim) = diff_centers(1:sb%dim) - shift(1:sb%dim)

    !As we take the union of the two submeshes, we know that we have all the points from the first one included.
    !The extra points from the second submesh are those which are not included in the first one
    !At the moment np_part extra points are not included
    is = sm1%np
    do ip = 1, sm2%np
      !sm2%x contains points coordinates defined with respect to sm2%center
      xx(1:sb%dim) = sm2%x(ip, 1:sb%dim)-diff_centers(1:sb%dim)
      !If the point is not in sm1, we add it
      if(sum(xx(1:sb%dim)**2) > sm1%radius**2) is = is + 1
    end do 

    this%np = is
    this%np_part = this%np

    SAFE_ALLOCATE(this%map(1:this%np_part))
    SAFE_ALLOCATE(this%x(1:this%np_part, 0:sb%dim))
    this%map(1:sm1%np) = sm1%map(1:sm1%np)
    this%x(1:sm1%np, 0:sb%dim) = sm1%x(1:sm1%np, 0:sb%dim)

    !iterate again to fill the tables
    is = sm1%np
    do ip = 1, sm2%np
      xx(1:sb%dim) = sm2%x(ip, 1:sb%dim) - diff_centers(1:sb%dim)
      r2 = sum(xx(1:sb%dim)**2)
      if(r2 > sm1%radius**2) then
        is = is + 1
        this%map(is) = sm2%map(ip)
        this%x(is, 0) = sqrt(r2)
        this%x(is, 1:sb%dim) = xx(1:sb%dim)
      end if
    end do

    call profiling_out(prof)
    POP_SUB(submesh_merge)
  end subroutine submesh_merge

  ! --------------------------------------------------------------
  !This routine shifts the center of a submesh, without changing the grid points
  subroutine submesh_shift_center(this, sb, newcenter)
    type(submesh_t),      intent(inout)  :: this 
    type(simul_box_t),    intent(in)     :: sb
    FLOAT,                intent(in)     :: newcenter(:)
    
    FLOAT :: r2
    integer :: ip
    type(profile_t), save :: prof
    FLOAT :: xx(1:MAX_DIM), diff_centers(1:MAX_DIM), oldcenter(1:MAX_DIM)
    
    PUSH_SUB(submesh_shift_center)
    call profiling_in(prof, "SUBMESH_SHIFT")

    oldcenter(1:sb%dim) = this%center(1:sb%dim)
    this%center(1:sb%dim)  = newcenter(1:sb%dim)
   

    diff_centers(1:sb%dim) = newcenter(1:sb%dim)-oldcenter(1:sb%dim)

    do ip = 1, this%np
      xx(1:sb%dim) = this%x(ip, 1:sb%dim) - diff_centers(1:sb%dim)
      r2 = sum(xx(1:sb%dim)**2)
      this%x(ip, 0) = sqrt(r2)
      this%x(ip, 1:sb%dim) = xx(1:sb%dim)
    end do

    call profiling_out(prof)
    POP_SUB(submesh_shift_center)
  end subroutine submesh_shift_center


  ! --------------------------------------------------------------

  subroutine submesh_broadcast(this, mesh, center, radius, root, mpi_grp)
    type(submesh_t),      intent(inout)  :: this
    type(mesh_t), target, intent(in)     :: mesh
    FLOAT,                intent(in)     :: center(:)
    FLOAT,                intent(in)     :: radius
    integer,              intent(in)     :: root
    type(mpi_grp_t),      intent(in)     :: mpi_grp

    integer :: nparray(1:3)
    type(profile_t), save :: prof

    PUSH_SUB(submesh_broadcast)
    call profiling_in(prof, 'SUBMESH_BCAST')
    
    if(root /= mpi_grp%rank) then    
      this%mesh => mesh
      this%center(1:mesh%sb%dim) = center(1:mesh%sb%dim)
      this%radius = radius
    end if

    if(mpi_grp%size > 1) then

      if(root == mpi_grp%rank) then
        nparray(1) = this%np
        nparray(2) = this%np_part
        if(this%overlap) then 
          nparray(3) = 1
        else
          nparray(3) = 0
        end if
      end if

#ifdef HAVE_MPI
      call MPI_Bcast(nparray, 3, MPI_INTEGER, root, mpi_grp%comm, mpi_err)
      call MPI_Barrier(mpi_grp%comm, mpi_err)
#endif
      this%np = nparray(1)
      this%np_part = nparray(2)
      this%overlap = (nparray(3) == 1)

      if(root /= mpi_grp%rank) then
        SAFE_ALLOCATE(this%map(1:this%np_part))
        SAFE_ALLOCATE(this%x(1:this%np_part, 0:mesh%sb%dim))
      end if

#ifdef HAVE_MPI
      if(this%np_part > 0) then
        call MPI_Bcast(this%map(1), this%np_part, MPI_INTEGER, root, mpi_grp%comm, mpi_err)
        call MPI_Barrier(mpi_grp%comm, mpi_err)
        call MPI_Bcast(this%x(1, 0), this%np_part*(mesh%sb%dim + 1), MPI_FLOAT, root, mpi_grp%comm, mpi_err)
        call MPI_Barrier(mpi_grp%comm, mpi_err)
      end if
#endif

    end if

    call profiling_out(prof)
    POP_SUB(submesh_broadcast)
  end subroutine submesh_broadcast
   
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

  logical function submesh_overlap(sm1, sm2) result(overlap)
    type(submesh_t),      intent(in)   :: sm1
    type(submesh_t),      intent(in)   :: sm2
    
    integer :: ii, jj, dd
    FLOAT :: distance

    !no PUSH_SUB, called too often

    if(.not. simul_box_is_periodic(sm1%mesh%sb)) then
      !first check the distance
      distance = sum((sm1%center(1:sm1%mesh%sb%dim) - sm2%center(1:sm2%mesh%sb%dim))**2)
      overlap = distance <= (CNST(1.5)*(sm1%radius + sm2%radius))**2
      
      ! if they are very far, no need to check in detail
      if(.not. overlap) return
    end if
    
    ! Otherwise check whether they have the some point in common. We
    ! can make the comparison faster using that the arrays are sorted.
    overlap = .false.
    ii = 1
    jj = 1
    do while(ii <= sm1%np_part .and. jj <= sm2%np_part)
      dd = sm1%map(ii) - sm2%map(jj)
      if(dd < 0) then
        ii = ii + 1
      else if(dd > 0) then
        jj = jj + 1
      else
        overlap = .true.
        exit
      end if
    end do

#ifdef HAVE_MPI
    if(sm1%mesh%parallel_in_domains) then
      call MPI_Allreduce(MPI_IN_PLACE, overlap, 1, MPI_LOGICAL, MPI_LOR, sm1%mesh%mpi_grp%comm, mpi_err)
    end if
#endif
    
  end function submesh_overlap

  ! -------------------------------------------------------------

    subroutine submesh_build_global(this)
    type(submesh_t),      intent(inout)   :: this

    integer, allocatable :: part_np(:)
    integer :: ipart, ind, ip

    PUSH_SUB(submesh_build_global)

    if(.not. this%mesh%parallel_in_domains) then
      POP_SUB(submesh_build_global)
      return
    end if 

    SAFE_ALLOCATE(part_np(this%mesh%vp%npart))
    part_np = 0
    part_np(this%mesh%vp%partno) = this%np

  #if defined(HAVE_MPI)
    call comm_allreduce(this%mesh%mpi_grp%comm, part_np)
  #endif 
    this%np_global = sum(part_np)

    SAFE_ALLOCATE(this%x_global(1:this%np_global, 1:this%mesh%sb%dim))
    SAFE_ALLOCATE(this%part_v(1:this%np_global))
    SAFE_ALLOCATE(this%global2local(1:this%np_global))
    this%x_global(1:this%np_global, 1:this%mesh%sb%dim) = M_ZERO
    this%part_v(1:this%np_global) = 0
    this%global2local(1:this%np_global) = 0

    ind = 0
    do ipart = 1, this%mesh%vp%npart
      if(ipart == this%mesh%vp%partno) then
        do ip = 1, this%np
          this%x_global(ind + ip, 1:this%mesh%sb%dim) = this%x(ip,1:this%mesh%sb%dim)
          this%part_v(ind + ip) = this%mesh%vp%partno
          this%global2local(ind + ip) = ip
        end do
      end if
      ind = ind + part_np(ipart)
    end do 

   #if defined(HAVE_MPI)
    call comm_allreduce(this%mesh%mpi_grp%comm, this%x_global)
    call comm_allreduce(this%mesh%mpi_grp%comm, this%part_v)
    call comm_allreduce(this%mesh%mpi_grp%comm, this%global2local)
   #endif 

    SAFE_DEALLOCATE_A(part_np)

    POP_SUB(submesh_build_global)
  end subroutine submesh_build_global


  subroutine submesh_end_global(this)
    type(submesh_t),      intent(inout)   :: this

    PUSH_SUB(submesh_end_global)

    SAFE_DEALLOCATE_A(this%x_global)
    this%np_global = -1
    SAFE_DEALLOCATE_A(this%part_v)
    SAFE_DEALLOCATE_A(this%global2local)

    POP_SUB(submesh_end_global)
  end subroutine submesh_end_global

  
  ! -----------------------------------------------------------
  
  subroutine zzsubmesh_add_to_mesh(this, sphi, phi, factor)
    type(submesh_t),  intent(in)    :: this
    CMPLX,            intent(in)    :: sphi(:)
    CMPLX,            intent(inout) :: phi(:)
    CMPLX,  optional, intent(in)    :: factor
    
    integer :: ip, m
    
    PUSH_SUB(zzdsubmesh_add_to_mesh)
   
    if(present(factor)) then
      !Loop unrolling inspired by BLAS axpy routine
      m = mod(this%np,4)
      do ip = 1, m
        phi(this%map(ip)) = phi(this%map(ip)) + factor*sphi(ip)
      end do
      if( this%np.ge.4) then
        do ip = m+1, this%np, 4
          phi(this%map(ip))   = phi(this%map(ip))   + factor*sphi(ip)
          phi(this%map(ip+1)) = phi(this%map(ip+1)) + factor*sphi(ip+1)
          phi(this%map(ip+2)) = phi(this%map(ip+2)) + factor*sphi(ip+2)
          phi(this%map(ip+3)) = phi(this%map(ip+3)) + factor*sphi(ip+3)
        end do
      end if
    else
      m = mod(this%np,4)
      do ip = 1, m
        phi(this%map(ip)) = phi(this%map(ip)) + sphi(ip)
      end do
      if( this%np.ge.4) then
        do ip = m+1, this%np, 4
          phi(this%map(ip))   = phi(this%map(ip))   + sphi(ip)
          phi(this%map(ip+1)) = phi(this%map(ip+1)) + sphi(ip+1)
          phi(this%map(ip+2)) = phi(this%map(ip+2)) + sphi(ip+2)
          phi(this%map(ip+3)) = phi(this%map(ip+3)) + sphi(ip+3)
        end do
      end if
    end if 
    
    POP_SUB(zzdsubmesh_add_to_mesh)
  end subroutine zzsubmesh_add_to_mesh


  !------------------------------------------------------------

  CMPLX function zzsubmesh_to_mesh_dotp(this, sphi, phi, reduce) result(dotp)
    type(submesh_t),   intent(in) :: this
    CMPLX,             intent(in) :: sphi(:)
    CMPLX,             intent(in) :: phi(:)
    logical, optional, intent(in) :: reduce
  
    integer :: is, m, ip
  
    PUSH_SUB(zzsubmesh_to_mesh_dotp)
  
    dotp = cmplx(M_ZERO, M_ZERO)
  
    if(this%mesh%use_curvilinear) then
      do is = 1, this%np
        dotp = dotp + this%mesh%vol_pp(this%map(is))*phi(this%map(is))*conjg(sphi(is))
      end do
    else
      m = mod(this%np,4)
      do ip = 1, m
        dotp = dotp + phi(this%map(ip))*conjg(sphi(ip))
      end do
      if( this%np.ge.4) then
        do ip = m+1, this%np, 4
          dotp = dotp + phi(this%map(ip))*conjg(sphi(ip))     &
                      + phi(this%map(ip+1))*conjg(sphi(ip+1)) &
                      + phi(this%map(ip+2))*conjg(sphi(ip+2)) &
                      + phi(this%map(ip+3))*conjg(sphi(ip+3))
        end do
      end if
        dotp = dotp*this%mesh%vol_pp(1)
    end if
  
    if(optional_default(reduce, .true.) .and. this%mesh%parallel_in_domains) then
      call profiling_in(C_PROFILING_SM_REDUCE, "SM_REDUCE")
      call comm_allreduce(this%mesh%vp%comm, dotp)
      call profiling_out(C_PROFILING_SM_REDUCE)
    end if 
 
    POP_SUB(zzsubmesh_to_mesh_dotp)
  end function zzsubmesh_to_mesh_dotp

  !------------------------------------------------------------
  !> finds the dimension of a box containing the submesh
  subroutine submesh_get_cube_dim(sm, db, dim)
    type(submesh_t),   intent(in)  :: sm
    integer,           intent(out) :: db(:)
    integer,           intent(in)  :: dim

    integer :: ip, idir
    FLOAT :: chi(1:MAX_DIM)

    PUSH_SUB(submesh_get_cube_dim)

    db = 1

    do ip = 1, sm%np
      !TODO: should be curvilinear_x2chi here instead
      chi(1:dim) = matmul(sm%x(ip,1:dim), sm%mesh%sb%klattice_primitive(1:dim, 1:dim))
      
      do idir = 1, dim
        db(idir) = max(db(idir), nint(abs(chi(idir))/sm%mesh%spacing(idir) + M_HALF))
      end do
    end do

    do idir = 1, dim
      db(idir) = 2 * db(idir) + 1
    end do

    POP_SUB(submesh_get_cube_dim)
  end subroutine submesh_get_cube_dim

  !------------------------------------------------------------
  subroutine submesh_init_cube_map(sm, db, dim)
    type(submesh_t),   intent(inout)  :: sm
    integer,           intent(in)     :: db(:)
    integer,           intent(in)     :: dim

    integer :: ip, idir
    FLOAT :: chi(1:MAX_DIM), shift(1:MAX_DIM)

    PUSH_SUB(submesh_init_cube_map)

    sm%cube_map%nmap = sm%np
    SAFE_ALLOCATE(sm%cube_map%map(1:dim, sm%cube_map%nmap))

    !The center of the submesh does not belong to the mesh
    !So we first need to find the closest grid point, and center the cube to it
    chi(1:dim) = matmul(sm%center(1:dim), sm%mesh%sb%klattice_primitive(1:dim, 1:dim))
    do idir = 1, dim
      shift(idir) = nint(chi(idir)/sm%mesh%spacing(idir))*sm%mesh%spacing(idir)
    end do
    shift(1:dim) = matmul(sm%mesh%sb%rlattice_primitive(1:dim,1:dim), shift(1:dim)) 
    shift(1:dim) = shift(1:dim) - sm%center(1:dim) 

    do ip = 1, sm%np
      !TODO: should be curvilinear_x2chi here instead
      chi(1:dim) = matmul(sm%x(ip,1:dim) - shift(1:dim), sm%mesh%sb%klattice_primitive(1:dim, 1:dim))
      do idir = 1, dim
        sm%cube_map%map(idir,ip) = nint(chi(idir)/sm%mesh%spacing(idir))
      end do
    end do

    if(accel_is_enabled()) then
      call accel_create_buffer(sm%cube_map%map_buffer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, sm%cube_map%nmap*5)
      call accel_write_buffer(sm%cube_map%map_buffer, sm%cube_map%nmap*5, sm%cube_map%map)
    end if


    POP_SUB(submesh_init_cube_map)
  end subroutine submesh_init_cube_map

  !------------------------------------------------------------
  subroutine submesh_end_cube_map(sm)
    type(submesh_t),   intent(inout)  :: sm

    PUSH_SUB(submesh_end_cube_map)

    call mesh_cube_map_end(sm%cube_map)

    POP_SUB(submesh_end_cube_map)
  end subroutine submesh_end_cube_map


#include "undef.F90"
#include "real.F90"
#include "submesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "submesh_inc.F90"

end module submesh_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
