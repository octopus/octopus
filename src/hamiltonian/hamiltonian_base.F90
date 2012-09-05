!! Copyright (C) 2009 X. Andrade
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
!! $Id: hamiltonian_base.F90 3988 2008-03-31 15:06:50Z fnog $

#include "global.h"

module hamiltonian_base_m
  use batch_m
  ! do not include blas, since we pass complex values to dgemm
  ! use blas_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use datasets_m
  use derivatives_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_m
  use kb_projector_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use opencl_m
  use parser_m
  use profiling_m
  use projector_m
  use projector_matrix_m
  use ps_m
  use simul_box_m
  use species_m
  use states_m
  use states_dim_m
  use submesh_m
  use types_m
  use varinfo_m

  implicit none

  private

  public ::                                    &
    hamiltonian_base_t,                        &
    dhamiltonian_base_local,                   &
    zhamiltonian_base_local,                   &
    dhamiltonian_base_magnetic,                &
    zhamiltonian_base_magnetic,                &
    dhamiltonian_base_nlocal_start,            &
    zhamiltonian_base_nlocal_start,            &
    dhamiltonian_base_nlocal_finish,           &
    zhamiltonian_base_nlocal_finish,           &
    hamiltonian_base_has_magnetic,             &
    hamiltonian_base_init,                     &
    hamiltonian_base_end,                      &
    hamiltonian_base_allocate,                 &
    hamiltonian_base_clear,                    &
    hamiltonian_base_build_proj,               &
    hamiltonian_base_update,                   &
    projection_t

  !> This object stores and applies an electromagnetic potential that
  !! can be represented by different types of potentials.

  type hamiltonian_base_t
    integer                           :: nspin
    type(nl_operator_t),      pointer :: kinetic
    type(projector_matrix_t), pointer :: projector_matrices(:) 
    FLOAT,                    pointer :: potential(:, :)
    FLOAT,                    pointer :: Impotential(:, :)!cmplxscl
    FLOAT,                    pointer :: uniform_magnetic_field(:)
    FLOAT,                    pointer :: uniform_vector_potential(:)
    FLOAT,                    pointer :: vector_potential(:, :)
    integer                           :: nprojector_matrices
    logical                           :: apply_projector_matrices
    integer                           :: full_projection_size
    integer                           :: max_npoints
    integer                           :: total_points
    integer                           :: max_nprojs
    CMPLX,                    pointer :: projector_phases(:, :, :)
    integer,                  pointer :: projector_to_atom(:)
    integer                           :: nregions
    integer,                  pointer :: regions(:)
#ifdef HAVE_OPENCL
    type(opencl_mem_t)                :: potential_opencl
    type(opencl_mem_t)                :: buff_offsets
    type(opencl_mem_t)                :: buff_matrices
    type(opencl_mem_t)                :: buff_maps
    type(opencl_mem_t)                :: buff_scals
    type(opencl_mem_t)                :: buff_pos
    type(opencl_mem_t)                :: buff_invmap
#endif
  end type hamiltonian_base_t

  type projection_t
    FLOAT, pointer     :: dprojection(:, :)
    CMPLX, pointer     :: zprojection(:, :)
#ifdef HAVE_OPENCL
    type(opencl_mem_t) :: buff_projection
#endif
  end type projection_t

  integer, public ::                     &
    TERM_ALL                 = HUGE(1),  &
    TERM_KINETIC             =   1,      &
    TERM_LOCAL_POTENTIAL     =   2,      & 
    TERM_NON_LOCAL_POTENTIAL =   4,      &
    TERM_OTHERS              =   8,      &
    TERM_LOCAL_EXTERNAL      =  16,      &
    TERM_MGGA                =  32

  integer, public ::                       &
    FIELD_POTENTIAL                = 1,    &
    FIELD_VECTOR_POTENTIAL         = 2,    &
    FIELD_UNIFORM_VECTOR_POTENTIAL = 4,    &
    FIELD_UNIFORM_MAGNETIC_FIELD   = 8
  

  type(profile_t), save :: prof_vnlpsi_start, prof_vnlpsi_finish, prof_magnetic, prof_vlpsi, prof_gather, prof_scatter

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_init(this, nspin)
    type(hamiltonian_base_t), intent(inout) :: this
    integer,                  intent(in)    :: nspin

    PUSH_SUB(hamiltonian_base_init)

    this%nspin = nspin

    nullify(this%potential)
    nullify(this%Impotential)!cmplxscl
    nullify(this%uniform_magnetic_field)
    nullify(this%uniform_vector_potential)
    nullify(this%vector_potential)
    nullify(this%projector_matrices)
    nullify(this%regions)
    this%apply_projector_matrices = .false.
    this%nprojector_matrices = 0

    POP_SUB(hamiltonian_base_init)
  end subroutine hamiltonian_base_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_end(this)
    type(hamiltonian_base_t), intent(inout) :: this

    PUSH_SUB(hamiltonian_base_end)

    SAFE_DEALLOCATE_P(this%potential)
    SAFE_DEALLOCATE_P(this%Impotential)!cmplxscl
    SAFE_DEALLOCATE_P(this%vector_potential)
    SAFE_DEALLOCATE_P(this%uniform_vector_potential)
    SAFE_DEALLOCATE_P(this%uniform_magnetic_field)
    call hamiltonian_base_destroy_proj(this)

    POP_SUB(hamiltonian_base_end)
  end subroutine hamiltonian_base_end

  ! ---------------------------------------------------------- 
  !
  !> This functions sets to zero all fields that are currently
  !! allocated.
  !
  subroutine hamiltonian_base_clear(this)
    type(hamiltonian_base_t), intent(inout) :: this

    PUSH_SUB(hamiltonian_clear)

    if(associated(this%potential))                this%potential = M_ZERO
    if(associated(this%Impotential))              this%Impotential = M_ZERO!cmplxscl
    if(associated(this%uniform_vector_potential)) this%uniform_vector_potential = M_ZERO
    if(associated(this%vector_potential))         this%vector_potential = M_ZERO
    if(associated(this%uniform_magnetic_field))   this%uniform_magnetic_field = M_ZERO

    POP_SUB(hamiltonian_clear)
  end subroutine hamiltonian_base_clear


  ! ---------------------------------------------------------------
  !> This function ensures that the corresponding field is allocated.
  subroutine hamiltonian_base_allocate(this, mesh, field, cmplxscl)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: field
    logical,                  intent(in)    :: cmplxscl

    PUSH_SUB(hamiltonian_base_allocate)

    if(iand(FIELD_POTENTIAL, field) /= 0) then 
      if(.not. associated(this%potential)) then
        SAFE_ALLOCATE(this%potential(1:mesh%np, 1:this%nspin))
        this%potential = M_ZERO
        if(cmplxscl) then
          SAFE_ALLOCATE(this%Impotential(1:mesh%np, 1:this%nspin))
          this%Impotential = M_ZERO
        end if
#ifdef HAVE_OPENCL
        if(opencl_is_enabled()) then
          call opencl_create_buffer(this%potential_opencl, CL_MEM_READ_ONLY, TYPE_FLOAT, opencl_padded_size(mesh%np)*this%nspin)
        end if
#endif
      end if
    end if

    if(iand(FIELD_UNIFORM_VECTOR_POTENTIAL, field) /= 0) then 
      if(.not. associated(this%uniform_vector_potential)) then
        SAFE_ALLOCATE(this%uniform_vector_potential(1:MAX_DIM))
        this%uniform_vector_potential = M_ZERO
      end if
    end if

    if(iand(FIELD_VECTOR_POTENTIAL, field) /= 0) then 
      if(.not. associated(this%vector_potential)) then
        SAFE_ALLOCATE(this%vector_potential(1:MAX_DIM, 1:mesh%np))
        this%vector_potential = M_ZERO
      end if
    end if

    if(iand(FIELD_UNIFORM_MAGNETIC_FIELD, field) /= 0) then 
      if(.not. associated(this%uniform_magnetic_field)) then
        SAFE_ALLOCATE(this%uniform_magnetic_field(1:MAX_DIM))
        this%uniform_magnetic_field = M_ZERO
      end if
    end if

    POP_SUB(hamiltonian_base_allocate)
  end subroutine hamiltonian_base_allocate

  ! ---------------------------------------------------------- 
  !
  !> If both a uniform and non-uniform vector potentials are allocated,
  !! this function copies the uniform in the non-uniform one. In the
  !! future it may perform other internal consistency operations.
  !
  subroutine hamiltonian_base_update(this, mesh)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh

    integer :: ispin
    integer :: offset

    PUSH_SUB(hamiltonian_base_update)

    if(associated(this%uniform_vector_potential) .and. associated(this%vector_potential)) then
      call unify_vector_potentials()
    end if

#ifdef HAVE_OPENCL
    if(associated(this%potential) .and. opencl_is_enabled()) then

      offset = 0
      do ispin = 1, this%nspin
        call opencl_write_buffer(this%potential_opencl, mesh%np, this%potential(:, ispin), offset = offset)
        offset = offset + opencl_padded_size(mesh%np)
      end do

    end if
#endif

    POP_SUB(hamiltonian_base_update)

  contains

    subroutine unify_vector_potentials()
      integer :: idir, ip

      PUSH_SUB(hamiltonian_base_update.unify_vector_potentials)
      
      ! copy the uniform vector potential onto the non-uniform one
      forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np) 
        this%vector_potential(idir, ip) = &
          this%vector_potential(idir, ip) + this%uniform_vector_potential(idir)
      end forall
      
      ! and deallocate
      SAFE_DEALLOCATE_P(this%uniform_vector_potential)
      POP_SUB(hamiltonian_base_update.unify_vector_potentials)      
    end subroutine unify_vector_potentials

  end subroutine hamiltonian_base_update
  
  !--------------------------------------------------------

  subroutine hamiltonian_base_destroy_proj(this)
    type(hamiltonian_base_t), intent(inout) :: this

    integer :: iproj

    PUSH_SUB(hamiltonian_base_destroy_proj)

    if(associated(this%projector_matrices)) then

#ifdef HAVE_OPENCL
      if(opencl_is_enabled()) then
        call opencl_release_buffer(this%buff_offsets)
        call opencl_release_buffer(this%buff_matrices)
        call opencl_release_buffer(this%buff_maps)
        call opencl_release_buffer(this%buff_scals)
        call opencl_release_buffer(this%buff_pos)
        call opencl_release_buffer(this%buff_invmap)
      end if
#endif

      do iproj = 1, this%nprojector_matrices
        call projector_matrix_deallocate(this%projector_matrices(iproj))
      end do
      SAFE_DEALLOCATE_P(this%regions)
      SAFE_DEALLOCATE_P(this%projector_matrices)
      SAFE_DEALLOCATE_P(this%projector_phases)
      SAFE_DEALLOCATE_P(this%projector_to_atom)
    end if

    POP_SUB(hamiltonian_base_destroy_proj)
  end subroutine hamiltonian_base_destroy_proj

  !-----------------------------------------------------------------
    
  subroutine hamiltonian_base_build_proj(this, mesh, epot, geo)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    type(epot_t),             intent(in)    :: epot
    type(geometry_t),         intent(in)    :: geo

    integer :: iatom, iproj, ll, lmax, lloc, mm, ic
    integer :: nmat, imat, ip, iorder
    integer :: nregion, jatom, katom, iregion
    integer, allocatable :: order(:), head(:), region_count(:)
    logical, allocatable :: atom_counted(:)
    logical :: overlap
    type(ps_t), pointer :: ps
    type(projector_matrix_t), pointer :: pmat
    type(kb_projector_t),     pointer :: kb_p
    type(profile_t), save :: color_prof

    PUSH_SUB(hamiltonian_base_build_proj)

    call profiling_in(color_prof, "ATOM_COLORING")

    ! this is most likely a very inefficient algorithm, O(natom**2) or
    ! O(natom**3), probably it should be replaced by something better.

    SAFE_ALLOCATE(order(1:epot%natoms))
    SAFE_ALLOCATE(head(1:epot%natoms + 1))
    SAFE_ALLOCATE(region_count(1:epot%natoms))
    SAFE_ALLOCATE(atom_counted(1:epot%natoms))

    atom_counted = .false.
    order = -1

    head(1) = 1
    nregion = 0
    do 
      nregion = nregion + 1
      ASSERT(nregion <= epot%natoms)

      region_count(nregion) = 0

      do iatom = 1, epot%natoms
        if(atom_counted(iatom)) cycle

        overlap = .false.

        if(projector_is(epot%proj(iatom), M_KB)) then
          do jatom = 1, region_count(nregion)
            katom = order(head(nregion) + jatom - 1)
            overlap = submesh_overlap(epot%proj(iatom)%sphere, epot%proj(katom)%sphere)
            if(overlap) exit
          end do
        end if

        if(.not. overlap) then
          INCR(region_count(nregion), 1)
          order(head(nregion) - 1 + region_count(nregion)) = iatom
          atom_counted(iatom) = .true.
        end if

      end do

      head(nregion + 1) = head(nregion) + region_count(nregion)

      if(all(atom_counted)) exit
    end do

    if(in_debug_mode) then
      call messages_write('The atoms can be separated in ')
      call messages_write(nregion)
      call messages_write(' non-overlapping groups.')
      call messages_info()
    end if

    do iregion = 1, nregion
      do iatom = head(iregion), head(iregion + 1) - 1
        if(.not. projector_is(epot%proj(order(iatom)), M_KB)) cycle
        do jatom = head(iregion), iatom - 1
          if(.not. projector_is(epot%proj(order(jatom)), M_KB)) cycle
          ASSERT(.not. submesh_overlap(epot%proj(order(iatom))%sphere, epot%proj(order(jatom))%sphere))
        end do
      end do
    end do

    call profiling_out(color_prof)

    ! deallocate previous projectors
    call hamiltonian_base_destroy_proj(this)

    ! count projectors
    this%nprojector_matrices = 0
    this%apply_projector_matrices = .false.
    this%nregions = nregion

    do iorder = 1, epot%natoms
      iatom = order(iorder)

      if(projector_is(epot%proj(iatom), M_KB)) then
        INCR(this%nprojector_matrices, 1)
        this%apply_projector_matrices = .true.
      else if(.not. projector_is_null(epot%proj(iatom))) then
        ! for the moment only KB projectors are supported
        this%apply_projector_matrices = .false.
        exit
      end if
    end do

    if(mesh%use_curvilinear) this%apply_projector_matrices = .false.
    if(simul_box_is_periodic(mesh%sb) .and. opencl_is_enabled()) this%apply_projector_matrices = .false.

    if(.not. this%apply_projector_matrices) then
      POP_SUB(hamiltonian_base_build_proj)
      return
    end if

    SAFE_ALLOCATE(this%projector_matrices(1:this%nprojector_matrices))
    SAFE_ALLOCATE(this%regions(1:this%nprojector_matrices + 1))
    SAFE_ALLOCATE(this%projector_to_atom(1:epot%natoms))

    this%full_projection_size = 0
    this%regions(this%nregions + 1) = this%nprojector_matrices + 1

    iproj = 0
    do iregion = 1, this%nregions
      this%regions(iregion) = iproj + 1
      do iorder = head(iregion), head(iregion + 1) - 1

        iatom = order(iorder)

        if(.not. projector_is(epot%proj(iatom), M_KB)) cycle
        INCR(iproj, 1)

        this%projector_to_atom(iproj) = iatom

        lmax = epot%proj(iatom)%lmax
        lloc = epot%proj(iatom)%lloc

        ! count the number of projectors for this matrix
        nmat = 0
        do ll = 0, lmax
          if (ll == lloc) cycle
          do mm = -ll, ll
            INCR(nmat, epot%proj(iatom)%kb_p(ll, mm)%n_c)
          end do
        end do

        ps => species_ps(geo%atom(iatom)%spec)
        pmat => this%projector_matrices(iproj)

        call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%np, nmat)

        ! generate the matrix
        pmat%projectors = M_ZERO

        imat = 1
        do ll = 0, lmax
          if (ll == lloc) cycle
          do mm = -ll, ll
            kb_p =>  epot%proj(iatom)%kb_p(ll, mm)
            do ic = 1, kb_p%n_c
              forall(ip = 1:pmat%npoints) pmat%projectors(ip, imat) = kb_p%p(ip, ic)
              pmat%scal(imat) = kb_p%e(ic)*mesh%vol_pp(1)
              INCR(imat, 1)
            end do
          end do
        end do

        forall(ip = 1:pmat%npoints) pmat%map(ip) = epot%proj(iatom)%sphere%map(ip)

        INCR(this%full_projection_size, pmat%nprojs)

      end do
    end do

!    do iregion = 1, this%nregions
!      print*, iregion, this%regions(iregion), this%regions(iregion + 1) - 1
!    end do

    this%total_points = 0
    this%max_npoints = 0
    this%max_nprojs = 0
    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)

      this%max_npoints = max(this%max_npoints, pmat%npoints)
      this%max_nprojs = max(this%max_nprojs, pmat%nprojs)
      INCR(this%total_points, pmat%npoints)
    end do

    nullify(this%projector_phases)

#ifdef HAVE_OPENCL
    if(opencl_is_enabled()) call build_opencl()
#endif

    POP_SUB(hamiltonian_base_build_proj)
  contains

    subroutine build_opencl()
#ifdef HAVE_OPENCL
      integer              :: matrix_size, scal_size
      integer, allocatable :: cnt(:), invmap(:, :), invmap2(:), pos(:)
      integer, allocatable :: offsets(:, :)
      integer, parameter   :: POINTS = 1, PROJS = 2, MATRIX = 3, MAP = 4, SCAL = 5
      integer              :: ip, is, ii, ipos

      PUSH_SUB(hamiltonian_base_build_proj.build_opencl)

      SAFE_ALLOCATE(offsets(1:5, 1:this%nprojector_matrices))
      SAFE_ALLOCATE(cnt(1:mesh%np))

      cnt = 0

      ! first we count
      matrix_size = 0
      this%total_points = 0
      scal_size = 0
      this%max_npoints = 0
      this%max_nprojs = 0
      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)

        this%max_npoints = max(this%max_npoints, pmat%npoints)
        this%max_nprojs = max(this%max_nprojs, pmat%nprojs)

        offsets(POINTS, imat) = pmat%npoints
        offsets(PROJS, imat) = pmat%nprojs

        offsets(MATRIX, imat) = matrix_size
        INCR(matrix_size, pmat%npoints*pmat%nprojs)

        offsets(MAP, imat) = this%total_points
        INCR(this%total_points, pmat%npoints)

        offsets(SCAL, imat) = scal_size
        INCR(scal_size, pmat%nprojs)

        do is = 1, pmat%npoints
          ip = pmat%map(is)
          INCR(cnt(ip), 1)
        end do
      end do

      SAFE_ALLOCATE(invmap(1:maxval(cnt), 1:mesh%np))
      SAFE_ALLOCATE(invmap2(1:maxval(cnt)*mesh%np))
      SAFE_ALLOCATE(pos(1:mesh%np + 1))

      cnt = 0
      ii = 0
      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)
        do is = 1, pmat%npoints
          ip = pmat%map(is)
          INCR(cnt(ip), 1)
          invmap(cnt(ip), ip) = ii
          INCR(ii, 1)
        end do
      end do

      ipos = 0
      pos(1) = 0
      do ip = 1, mesh%np
        do ii = 1, cnt(ip)
          INCR(ipos, 1)
          invmap2(ipos) = invmap(ii, ip)
        end do
        pos(ip + 1) = ipos
      end do

      ! allocate
      call opencl_create_buffer(this%buff_matrices, CL_MEM_READ_ONLY, TYPE_FLOAT, matrix_size)
      call opencl_create_buffer(this%buff_maps, CL_MEM_READ_ONLY, TYPE_INTEGER, this%total_points)
      call opencl_create_buffer(this%buff_scals, CL_MEM_READ_ONLY, TYPE_FLOAT, scal_size)

      ! now copy
      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)

        call opencl_write_buffer(this%buff_matrices, pmat%nprojs*pmat%npoints, pmat%projectors, offset = offsets(MATRIX, imat))
        call opencl_write_buffer(this%buff_maps, pmat%npoints, pmat%map, offset = offsets(MAP, imat))
        call opencl_write_buffer(this%buff_scals, pmat%nprojs, pmat%scal, offset = offsets(SCAL, imat))
      end do

      ! write the offsets
      call opencl_create_buffer(this%buff_offsets, CL_MEM_READ_ONLY, TYPE_INTEGER, 5*this%nprojector_matrices)
      call opencl_write_buffer(this%buff_offsets, 5*this%nprojector_matrices, offsets)

      ! the inverse map
      call opencl_create_buffer(this%buff_pos, CL_MEM_READ_ONLY, TYPE_INTEGER, mesh%np + 1)
      call opencl_write_buffer(this%buff_pos, mesh%np + 1, pos)

      call opencl_create_buffer(this%buff_invmap, CL_MEM_READ_ONLY, TYPE_INTEGER, ipos)
      call opencl_write_buffer(this%buff_invmap, ipos, invmap2)

      SAFE_DEALLOCATE_A(offsets)
      SAFE_DEALLOCATE_A(cnt)
      SAFE_DEALLOCATE_A(invmap)
      SAFE_DEALLOCATE_A(invmap2)
      SAFE_DEALLOCATE_A(pos)

#endif  

      POP_SUB(hamiltonian_base_build_proj.build_opencl)
    end subroutine build_opencl

  end subroutine hamiltonian_base_build_proj
    
  ! ----------------------------------------------------------------------------------

  logical pure function hamiltonian_base_has_magnetic(this) result(has_magnetic)
    type(hamiltonian_base_t), intent(in) :: this
    
    has_magnetic = associated(this%vector_potential) &
      .or. associated(this%uniform_magnetic_field)
    
  end function hamiltonian_base_has_magnetic

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_base_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_base_inc.F90"

end module hamiltonian_base_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
