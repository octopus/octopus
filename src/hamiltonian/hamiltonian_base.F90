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
  use blas_m
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
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use nl_operator_m
#ifdef HAVE_OPENCL
  use opencl_m
#endif
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
    dhamiltonian_base_non_local,               &
    zhamiltonian_base_non_local,               &
    hamiltonian_base_has_magnetic,             &
    hamiltonian_base_init,                     &
    hamiltonian_base_end,                      &
    hamiltonian_base_allocate,                 &
    hamiltonian_base_clear,                    &
    hamiltonian_base_build_proj,               &
    hamiltonian_base_update

  ! This object stores and applies an electromagnetic potential that
  ! can be represented by different types of potentials.

  type hamiltonian_base_t
    integer                           :: nspin
    type(nl_operator_t),      pointer :: kinetic
    type(projector_matrix_t), pointer :: projector_matrices(:) 
    FLOAT,                    pointer :: potential(:, :)
    FLOAT,                    pointer :: uniform_magnetic_field(:)
    FLOAT,                    pointer :: uniform_vector_potential(:)
    FLOAT,                    pointer :: vector_potential(:, :)
    integer                           :: nprojector_matrices
    logical                           :: apply_projector_matrices
#ifdef HAVE_OPENCL
    type(opencl_mem_t)           :: potential_opencl
#endif
  end type hamiltonian_base_t

  integer, public ::                     &
    TERM_ALL                 = HUGE(1),  &
    TERM_KINETIC             =   1,      &
    TERM_LOCAL_POTENTIAL     =   2,      & 
    TERM_NON_LOCAL_POTENTIAL =   4,      &
    TERM_OTHERS              =   8,      &
    TERM_LOCAL_EXTERNAL      =  16

  integer, public ::                       &
    FIELD_POTENTIAL                = 1,    &
    FIELD_VECTOR_POTENTIAL         = 2,    &
    FIELD_UNIFORM_VECTOR_POTENTIAL = 4,    &
    FIELD_UNIFORM_MAGNETIC_FIELD   = 8
  

  type(profile_t), save :: prof_vlpsi, prof_magnetic, prof_vnlpsi, prof_gather, prof_scatter

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_init(this, mesh, nspin)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: nspin

    call push_sub('hamiltonian_base.hamiltonian_base_init')

    this%nspin = nspin

    nullify(this%potential)
    nullify(this%uniform_magnetic_field)
    nullify(this%uniform_vector_potential)
    nullify(this%vector_potential)
    nullify(this%projector_matrices)
    this%apply_projector_matrices = .false.
    this%nprojector_matrices = 0

    call pop_sub('hamiltonian_base.hamiltonian_base_init')
  end subroutine hamiltonian_base_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_end(this)
    type(hamiltonian_base_t), intent(inout) :: this

    call push_sub('hamiltonian_base.hamiltonian_base_end')

    SAFE_DEALLOCATE_P(this%potential)
    SAFE_DEALLOCATE_P(this%vector_potential)
    SAFE_DEALLOCATE_P(this%uniform_vector_potential)
    SAFE_DEALLOCATE_P(this%uniform_magnetic_field)
    call hamiltonian_base_destroy_proj(this)

    call pop_sub('hamiltonian_base.hamiltonian_base_end')
  end subroutine hamiltonian_base_end

  ! ---------------------------------------------------------- 
  !
  ! This functions sets to zero all fields that are currently
  ! allocated.
  !
  subroutine hamiltonian_base_clear(this)
    type(hamiltonian_base_t), intent(inout) :: this

    call push_sub('hamiltonian_base.hamiltonian_clear')

    if(associated(this%potential))                this%potential = M_ZERO
    if(associated(this%uniform_vector_potential)) this%uniform_vector_potential = M_ZERO
    if(associated(this%vector_potential))         this%vector_potential = M_ZERO
    if(associated(this%uniform_magnetic_field))   this%uniform_magnetic_field = M_ZERO

    call pop_sub('hamiltonian_base.hamiltonian_clear')
  end subroutine hamiltonian_base_clear


  ! ---------------------------------------------------------------
  ! This function ensures that the corresponding field is allocated.
  subroutine hamiltonian_base_allocate(this, mesh, field)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: field
    call push_sub('hamiltonian_base.hamiltonian_base_allocate')

    if(iand(FIELD_POTENTIAL, field) /= 0) then 
      if(.not. associated(this%potential)) then
        SAFE_ALLOCATE(this%potential(1:mesh%np, 1:this%nspin))
        this%potential = M_ZERO
#ifdef HAVE_OPENCL
        if(opencl_is_available()) then
          call opencl_create_buffer(this%potential_opencl, CL_MEM_READ_ONLY, TYPE_FLOAT, opencl_padded_size(mesh%np))
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

    call pop_sub('hamiltonian_base.hamiltonian_base_allocate')
  end subroutine hamiltonian_base_allocate

  ! ---------------------------------------------------------- 
  !
  ! If both a uniform and non-uniform vector potentials are allocated,
  ! this function copies the uniform in the non-uniform one. In the
  ! future it may perform other internal consistency operations.
  !
  subroutine hamiltonian_base_update(this, mesh)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh

    call push_sub('hamiltonian_base.hamiltonian_check')

    if(associated(this%uniform_vector_potential) .and. associated(this%vector_potential)) then
      call unify_vector_potentials()
    end if

#ifdef HAVE_OPENCL
    if(associated(this%potential) .and. opencl_is_available()) then
      call opencl_write_buffer(this%potential_opencl, mesh%np, this%potential(:, 1))
    end if
#endif

    call pop_sub('hamiltonian_base.hamiltonian_check')

  contains

    subroutine unify_vector_potentials()
      integer :: idir, ip
      
      ! copy the uniform vector potential onto the non-uniform one
      forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np) 
        this%vector_potential(idir, ip) = &
          this%vector_potential(idir, ip) + this%uniform_vector_potential(idir)
      end forall
      
      ! and deallocate
      SAFE_DEALLOCATE_P(this%uniform_vector_potential)
      nullify(this%uniform_vector_potential)
      
    end subroutine unify_vector_potentials

  end subroutine hamiltonian_base_update
  
  !--------------------------------------------------------

  subroutine hamiltonian_base_destroy_proj(this)
    type(hamiltonian_base_t), intent(inout) :: this

    integer :: iproj
    
    if(associated(this%projector_matrices)) then
      do iproj = 1, this%nprojector_matrices
        call projector_matrix_deallocate(this%projector_matrices(1))
      end do
      SAFE_DEALLOCATE_P(this%projector_matrices)
    end if
  end subroutine hamiltonian_base_destroy_proj

  !-----------------------------------------------------------------
    
  subroutine hamiltonian_base_build_proj(this, mesh, epot, geo)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    type(epot_t),             intent(in)    :: epot
    type(geometry_t),         intent(in)    :: geo

    integer :: iatom, iproj, ll, lmax, lloc, mm, ic
    integer :: nmat, imat, ip
    type(ps_t), pointer :: ps
    type(projector_matrix_t), pointer :: pmat
    type(kb_projector_t),     pointer :: kb_p

    ! deallocate previous projectors
    call hamiltonian_base_destroy_proj(this)

    ! count projectors
    this%nprojector_matrices = 0
    this%apply_projector_matrices = .false.

    do iatom = 1, epot%natoms
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
    if(simul_box_is_periodic(mesh%sb)) this%apply_projector_matrices = .false.

    if(.not. this%apply_projector_matrices) return
    
    SAFE_ALLOCATE(this%projector_matrices(1:this%nprojector_matrices))

    iproj = 0
    do iatom = 1, epot%natoms
      if(.not. projector_is(epot%proj(iatom), M_KB)) cycle
      INCR(iproj, 1)

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

      call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%ns, nmat)

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

      forall(ip = 1:pmat%npoints) pmat%map(ip) = epot%proj(iatom)%sphere%jxyz(ip)

    end do
      
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
