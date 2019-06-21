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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module hamiltonian_base_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use geometry_oct_m
  use global_oct_m
  use hardware_oct_m
  use hgh_projector_oct_m
  use kb_projector_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use nl_operator_oct_m
  use profiling_oct_m
  use projector_oct_m
  use projector_matrix_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use types_oct_m

  implicit none

  private

  public ::                                    &
    hamiltonian_base_t,                        &
    dhamiltonian_base_local,                   &
    zhamiltonian_base_local,                   &
    dhamiltonian_base_local_sub,               &
    zhamiltonian_base_local_sub,               &
    dhamiltonian_base_magnetic,                &
    zhamiltonian_base_magnetic,                &
    dhamiltonian_base_rashba,                  &
    zhamiltonian_base_rashba,                  &
    dhamiltonian_base_nlocal_start,            &
    zhamiltonian_base_nlocal_start,            &
    dhamiltonian_base_nlocal_finish,           &
    zhamiltonian_base_nlocal_finish,           &
    dhamiltonian_base_nlocal_position_commutator, &
    zhamiltonian_base_nlocal_position_commutator, &
    hamiltonian_base_has_magnetic,             &
    hamiltonian_base_init,                     &
    hamiltonian_base_end,                      &
    hamiltonian_base_allocate,                 &
    hamiltonian_base_clear,                    &
    hamiltonian_base_build_proj,               &
    hamiltonian_base_update,                   &
    dhamiltonian_base_phase,                   &
    zhamiltonian_base_phase,                   &
    dhamiltonian_base_nlocal_force,            &
    zhamiltonian_base_nlocal_force,            &
    projection_t,                              &
    hamiltonian_base_projector_self_overlap

  !> This object stores and applies an electromagnetic potential that
  !! can be represented by different types of potentials.

  type hamiltonian_base_t
    private
    integer                                       :: nspin
    FLOAT                                         :: mass  !< Needed to compute the magnetic terms, if the mass is not one.
    FLOAT                                         :: rashba_coupling
    type(nl_operator_t),      pointer,     public :: kinetic
    type(projector_matrix_t), allocatable, public :: projector_matrices(:) 
    FLOAT,                    allocatable, public :: potential(:, :)
    FLOAT,                    allocatable, public :: Impotential(:, :)
    FLOAT,                    allocatable, public :: uniform_magnetic_field(:)
    FLOAT,                    allocatable, public :: uniform_vector_potential(:)
    FLOAT,                    allocatable, public :: vector_potential(:, :)
    integer,                               public :: nprojector_matrices
    logical,                               public :: apply_projector_matrices
    integer                                       :: full_projection_size
    integer,                               public :: max_npoints
    integer,                               public :: total_points
    integer                                       :: max_nprojs
    logical                                       :: projector_mix
    CMPLX,                    allocatable, public :: projector_phases(:, :, :)
    integer,                  allocatable, public :: projector_to_atom(:)
    integer                                       :: nregions
    integer,                  allocatable         :: regions(:)
    type(accel_mem_t)                             :: potential_opencl
    type(accel_mem_t)                             :: buff_offsets
    type(accel_mem_t)                             :: buff_matrices
    type(accel_mem_t)                             :: buff_maps
    type(accel_mem_t)                             :: buff_scals
    type(accel_mem_t)                             :: buff_position
    type(accel_mem_t)                             :: buff_pos
    type(accel_mem_t)                             :: buff_invmap
    type(accel_mem_t),                     public :: buff_projector_phases
    type(accel_mem_t)                             :: buff_mix
    CMPLX,                    pointer,     public :: phase(:, :)
    CMPLX,                    allocatable, public :: phase_corr(:,:)
    type(accel_mem_t),                     public :: buff_phase
    integer,                               public :: buff_phase_qn_start
    logical                                       :: projector_self_overlap  !< if .true. some projectors overlap with themselves
  end type hamiltonian_base_t

  type projection_t
    private
    FLOAT, allocatable     :: dprojection(:, :)
    CMPLX, allocatable     :: zprojection(:, :)
    type(accel_mem_t)      :: buff_projection
  end type projection_t

  integer, parameter, public ::          &
    TERM_ALL                 = HUGE(1),  &
    TERM_KINETIC             =   1,      &
    TERM_LOCAL_POTENTIAL     =   2,      & 
    TERM_NON_LOCAL_POTENTIAL =   4,      &
    TERM_OTHERS              =   8,      &
    TERM_LOCAL_EXTERNAL      =  16,      &
    TERM_MGGA                =  32,      &
    TERM_DFT_U               =  64 

  integer, parameter, public ::            &
    FIELD_POTENTIAL                = 1,    &
    FIELD_VECTOR_POTENTIAL         = 2,    &
    FIELD_UNIFORM_VECTOR_POTENTIAL = 4,    &
    FIELD_UNIFORM_MAGNETIC_FIELD   = 8

  type(profile_t), save :: prof_vnlpsi_start, prof_vnlpsi_finish, prof_magnetic, prof_vlpsi, prof_gather, prof_scatter, &
    prof_matelement, prof_matelement_gather, prof_matelement_reduce

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_init(this, nspin, mass, rashba_coupling)
    type(hamiltonian_base_t), intent(inout) :: this
    integer,                  intent(in)    :: nspin
    FLOAT,                    intent(in)    :: mass
    FLOAT,                    intent(in)    :: rashba_coupling

    PUSH_SUB(hamiltonian_base_init)

    this%nspin = nspin
    this%mass  = mass
    this%rashba_coupling = rashba_coupling

    this%apply_projector_matrices = .false.
    this%nprojector_matrices = 0

    this%projector_self_overlap = .false.
    
    POP_SUB(hamiltonian_base_init)
  end subroutine hamiltonian_base_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_end(this)
    type(hamiltonian_base_t), intent(inout) :: this

    PUSH_SUB(hamiltonian_base_end)

    if(allocated(this%potential) .and. accel_is_enabled()) then
      call accel_release_buffer(this%potential_opencl)
    end if
    
    SAFE_DEALLOCATE_A(this%potential)
    SAFE_DEALLOCATE_A(this%Impotential)
    SAFE_DEALLOCATE_A(this%vector_potential)
    SAFE_DEALLOCATE_A(this%uniform_vector_potential)
    SAFE_DEALLOCATE_A(this%uniform_magnetic_field)
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

    if(allocated(this%potential))                this%potential = M_ZERO
    if(allocated(this%Impotential))              this%Impotential = M_ZERO
    if(allocated(this%uniform_vector_potential)) this%uniform_vector_potential = M_ZERO
    if(allocated(this%vector_potential))         this%vector_potential = M_ZERO
    if(allocated(this%uniform_magnetic_field))   this%uniform_magnetic_field = M_ZERO

    POP_SUB(hamiltonian_clear)
  end subroutine hamiltonian_base_clear


  ! ---------------------------------------------------------------
  !> This function ensures that the corresponding field is allocated.
  subroutine hamiltonian_base_allocate(this, mesh, field, complex_potential)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: field
    logical,                  intent(in)    :: complex_potential

    PUSH_SUB(hamiltonian_base_allocate)

    if(bitand(FIELD_POTENTIAL, field) /= 0) then 
      if(.not. allocated(this%potential)) then
        SAFE_ALLOCATE(this%potential(1:mesh%np, 1:this%nspin))
        this%potential = M_ZERO
        if(complex_potential) then
          SAFE_ALLOCATE(this%Impotential(1:mesh%np, 1:this%nspin))
          this%Impotential = M_ZERO
        end if
        if(accel_is_enabled()) then
          call accel_create_buffer(this%potential_opencl, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, accel_padded_size(mesh%np)*this%nspin)
        end if
      end if
    end if

    if(bitand(FIELD_UNIFORM_VECTOR_POTENTIAL, field) /= 0) then 
      if(.not. allocated(this%uniform_vector_potential)) then
        SAFE_ALLOCATE(this%uniform_vector_potential(1:mesh%sb%dim))
        this%uniform_vector_potential = M_ZERO
      end if
    end if

    if(bitand(FIELD_VECTOR_POTENTIAL, field) /= 0) then 
      if(.not. allocated(this%vector_potential)) then
        SAFE_ALLOCATE(this%vector_potential(1:mesh%sb%dim, 1:mesh%np))
        this%vector_potential = M_ZERO
      end if
    end if

    if(bitand(FIELD_UNIFORM_MAGNETIC_FIELD, field) /= 0) then 
      if(.not. allocated(this%uniform_magnetic_field)) then
        SAFE_ALLOCATE(this%uniform_magnetic_field(1:max(mesh%sb%dim, 3)))
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

    if(allocated(this%uniform_vector_potential) .and. allocated(this%vector_potential)) then
      call unify_vector_potentials()
    end if

    if(allocated(this%potential) .and. accel_is_enabled()) then

      offset = 0
      do ispin = 1, this%nspin
        call accel_write_buffer(this%potential_opencl, mesh%np, this%potential(:, ispin), offset = offset)
        offset = offset + accel_padded_size(mesh%np)
      end do

    end if

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
      SAFE_DEALLOCATE_A(this%uniform_vector_potential)
      POP_SUB(hamiltonian_base_update.unify_vector_potentials)      
    end subroutine unify_vector_potentials

  end subroutine hamiltonian_base_update
  
  !--------------------------------------------------------

  subroutine hamiltonian_base_destroy_proj(this)
    type(hamiltonian_base_t), intent(inout) :: this

    integer :: iproj

    PUSH_SUB(hamiltonian_base_destroy_proj)

    if(allocated(this%projector_matrices)) then

      if(accel_is_enabled()) then
        call accel_release_buffer(this%buff_offsets)
        call accel_release_buffer(this%buff_matrices)
        call accel_release_buffer(this%buff_maps)
        call accel_release_buffer(this%buff_scals)
        call accel_release_buffer(this%buff_position)
        call accel_release_buffer(this%buff_pos)
        call accel_release_buffer(this%buff_invmap)
        if(this%projector_mix) call accel_release_buffer(this%buff_mix)
        if(allocated(this%projector_phases)) call accel_release_buffer(this%buff_projector_phases)
      end if

      do iproj = 1, this%nprojector_matrices
        call projector_matrix_deallocate(this%projector_matrices(iproj))
      end do
      SAFE_DEALLOCATE_A(this%regions)
      SAFE_DEALLOCATE_A(this%projector_matrices)
      SAFE_DEALLOCATE_A(this%projector_phases)
      SAFE_DEALLOCATE_A(this%projector_to_atom)
    end if

    POP_SUB(hamiltonian_base_destroy_proj)
  end subroutine hamiltonian_base_destroy_proj

  !-----------------------------------------------------------------
    
  subroutine hamiltonian_base_build_proj(this, mesh, epot)
    type(hamiltonian_base_t), target, intent(inout) :: this
    type(mesh_t),                     intent(in)    :: mesh
    type(epot_t),             target, intent(in)    :: epot

    integer :: iatom, iproj, ll, lmax, lloc, mm, ic, jc
    integer :: nmat, imat, ip, iorder
    integer :: nregion, jatom, katom, iregion
    integer, allocatable :: order(:), head(:), region_count(:)
    logical, allocatable :: atom_counted(:)
    logical :: overlap
    type(projector_matrix_t), pointer :: pmat
    type(kb_projector_t),     pointer :: kb_p
    type(hgh_projector_t),    pointer :: hgh_p
    type(profile_t), save :: color_prof

    PUSH_SUB(hamiltonian_base_build_proj)

    call profiling_in(color_prof, "ATOM_COLORING")

    ! this is most likely a very inefficient algorithm, O(natom**2) or
    ! O(natom**3), probably it should be replaced by something better.

    SAFE_ALLOCATE(order(1:epot%natoms))
    SAFE_ALLOCATE(head(1:epot%natoms + 1))
    SAFE_ALLOCATE(region_count(1:epot%natoms))
    SAFE_ALLOCATE(atom_counted(1:epot%natoms))

    this%projector_self_overlap = .false.
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

        if(.not. projector_is(epot%proj(iatom), PROJ_NONE)) then
          ASSERT(associated(epot%proj(iatom)%sphere%mesh))
          do jatom = 1, region_count(nregion)
            katom = order(head(nregion) + jatom - 1)
            if(projector_is(epot%proj(katom), PROJ_NONE)) cycle
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

    SAFE_DEALLOCATE_A(atom_counted)
    SAFE_DEALLOCATE_A(region_count)

    if(debug%info) then
      call messages_write('The atoms can be separated in ')
      call messages_write(nregion)
      call messages_write(' non-overlapping groups.')
      call messages_info()
    end if

    do iregion = 1, nregion
      do iatom = head(iregion), head(iregion + 1) - 1
        if(.not. projector_is(epot%proj(order(iatom)), PROJ_KB)) cycle
        do jatom = head(iregion), iatom - 1
          if(.not. projector_is(epot%proj(order(jatom)), PROJ_KB)) cycle
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

      if(projector_is(epot%proj(iatom), PROJ_KB) .or. projector_is(epot%proj(iatom), PROJ_HGH)) then
        INCR(this%nprojector_matrices, 1)
        this%apply_projector_matrices = .true.
      else if(.not. projector_is_null(epot%proj(iatom))) then
        this%apply_projector_matrices = .false.
        exit
      end if
    end do

    if(epot%reltype /= NOREL) this%apply_projector_matrices = .false.
    if(mesh%use_curvilinear)  this%apply_projector_matrices = .false.

    if(.not. this%apply_projector_matrices) then
      SAFE_DEALLOCATE_A(order)
      SAFE_DEALLOCATE_A(head)

      POP_SUB(hamiltonian_base_build_proj)
      return
    end if


    SAFE_ALLOCATE(this%projector_matrices(1:this%nprojector_matrices))
    SAFE_ALLOCATE(this%regions(1:this%nprojector_matrices + 1))
    SAFE_ALLOCATE(this%projector_to_atom(1:epot%natoms))

    this%full_projection_size = 0
    this%regions(this%nregions + 1) = this%nprojector_matrices + 1

    this%projector_mix = .false.
    
    iproj = 0
    do iregion = 1, this%nregions
      this%regions(iregion) = iproj + 1
      do iorder = head(iregion), head(iregion + 1) - 1

        iatom = order(iorder)

        if(projector_is(epot%proj(iatom), PROJ_NONE)) cycle
          
        INCR(iproj, 1)

        pmat => this%projector_matrices(iproj)

        this%projector_to_atom(iproj) = iatom

        lmax = epot%proj(iatom)%lmax
        lloc = epot%proj(iatom)%lloc

        if(projector_is(epot%proj(iatom), PROJ_KB)) then
          
          ! count the number of projectors for this matrix
          nmat = 0
          do ll = 0, lmax
            if (ll == lloc) cycle
            do mm = -ll, ll
              INCR(nmat, epot%proj(iatom)%kb_p(ll, mm)%n_c)
            end do
          end do
          
          call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%np, nmat, has_mix_matrix = .false.)
          
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

          this%projector_self_overlap = this%projector_self_overlap .or. epot%proj(iatom)%sphere%overlap

        else if(projector_is(epot%proj(iatom), PROJ_HGH)) then

          this%projector_mix = .true.
          
          ! count the number of projectors for this matrix
          nmat = 0
          do ll = 0, lmax
            if (ll == lloc) cycle
            do mm = -ll, ll
              nmat = nmat + 3
            end do
          end do
          
          call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%np, nmat, has_mix_matrix = .true.)

          ! generate the matrix
          pmat%projectors = M_ZERO
          pmat%mix = M_ZERO
          
          imat = 1
          do ll = 0, lmax
            if (ll == lloc) cycle
            do mm = -ll, ll
              hgh_p =>  epot%proj(iatom)%hgh_p(ll, mm)

              ! HGH pseudos mix different components, so we need to
              ! generate a matrix that mixes the projections
              do ic = 1, 3
                do jc = 1, 3
                  pmat%mix(imat - 1 + ic, imat - 1 + jc) = hgh_p%h(ic, jc)
                end do
              end do
              
              do ic = 1, 3
                forall(ip = 1:pmat%npoints) pmat%projectors(ip, imat) = hgh_p%dp(ip, ic)
                pmat%scal(imat) = mesh%volume_element
                INCR(imat, 1)
              end do
              
            end do
          end do

          this%projector_self_overlap = this%projector_self_overlap .or. epot%proj(iatom)%sphere%overlap
          
        else
          cycle          
        end if

        forall(ip = 1:pmat%npoints)
          pmat%map(ip) = epot%proj(iatom)%sphere%map(ip)
          pmat%position(1:3, ip) = epot%proj(iatom)%sphere%x(ip, 1:3)
        end forall

        INCR(this%full_projection_size, pmat%nprojs)

      end do
    end do

#ifdef HAVE_MPI
    if(mesh%parallel_in_domains) then
      call MPI_Allreduce(MPI_IN_PLACE, this%projector_self_overlap, 1, MPI_LOGICAL, MPI_LOR, mesh%mpi_grp%comm, mpi_err)
    end if
#endif
    
    SAFE_DEALLOCATE_A(order)
    SAFE_DEALLOCATE_A(head)

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

    if(accel_is_enabled()) call build_opencl()

    POP_SUB(hamiltonian_base_build_proj)

  contains

    subroutine build_opencl()
      integer              :: matrix_size, scal_size
      integer, allocatable :: cnt(:), invmap(:, :), invmap2(:), pos(:)
      integer, allocatable :: offsets(:, :)
      integer, parameter   :: OFFSET_SIZE = 6 ! also defined in share/opencl/projectors.cl
      integer, parameter   :: POINTS = 1, PROJS = 2, MATRIX = 3, MAP = 4, SCAL = 5, MIX = 6 ! update OFFSET_SIZE
      integer              :: ip, is, ii, ipos, mix_offset

      PUSH_SUB(hamiltonian_base_build_proj.build_opencl)

      SAFE_ALLOCATE(offsets(1:OFFSET_SIZE, 1:this%nprojector_matrices))
      SAFE_ALLOCATE(cnt(1:mesh%np))

      cnt = 0

      ! first we count
      matrix_size = 0
      this%total_points = 0
      scal_size = 0
      this%max_npoints = 0
      this%max_nprojs = 0
      mix_offset = 0
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

        if(allocated(pmat%mix)) then
          offsets(MIX, imat) = mix_offset
          INCR(mix_offset, pmat%nprojs**2)
        else
          offsets(MIX, imat) = -1
        end if
        
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
      call accel_create_buffer(this%buff_matrices, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, matrix_size)
      call accel_create_buffer(this%buff_maps, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%total_points)
      call accel_create_buffer(this%buff_position, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 3*this%total_points)
      call accel_create_buffer(this%buff_scals, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, scal_size)
      if(mix_offset > 0) call accel_create_buffer(this%buff_mix, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, mix_offset)
      
      ! now copy
      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)
        ! in parallel some spheres might not have points
        if(pmat%npoints > 0) then
          call accel_write_buffer(this%buff_matrices, pmat%nprojs*pmat%npoints, pmat%projectors, offset = offsets(MATRIX, imat))
          call accel_write_buffer(this%buff_maps, pmat%npoints, pmat%map, offset = offsets(MAP, imat))
          call accel_write_buffer(this%buff_position, 3*pmat%npoints, pmat%position, offset = 3*offsets(MAP, imat))
        end if
        call accel_write_buffer(this%buff_scals, pmat%nprojs, pmat%scal, offset = offsets(SCAL, imat))
        if(offsets(MIX, imat) /= -1) call accel_write_buffer(this%buff_mix, pmat%nprojs**2, pmat%mix, offset = offsets(MIX, imat))
      end do

      ! write the offsets
      call accel_create_buffer(this%buff_offsets, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, OFFSET_SIZE*this%nprojector_matrices)
      call accel_write_buffer(this%buff_offsets, OFFSET_SIZE*this%nprojector_matrices, offsets)

      ! the inverse map
      call accel_create_buffer(this%buff_pos, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, mesh%np + 1)
      call accel_write_buffer(this%buff_pos, mesh%np + 1, pos)

      call accel_create_buffer(this%buff_invmap, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, ipos)
      call accel_write_buffer(this%buff_invmap, ipos, invmap2)

      SAFE_DEALLOCATE_A(offsets)
      SAFE_DEALLOCATE_A(cnt)
      SAFE_DEALLOCATE_A(invmap)
      SAFE_DEALLOCATE_A(invmap2)
      SAFE_DEALLOCATE_A(pos)

      POP_SUB(hamiltonian_base_build_proj.build_opencl)
    end subroutine build_opencl

  end subroutine hamiltonian_base_build_proj
    
  ! ----------------------------------------------------------------------------------

  logical pure function hamiltonian_base_has_magnetic(this) result(has_magnetic)
    type(hamiltonian_base_t), intent(in) :: this
    
    has_magnetic = allocated(this%vector_potential) &
      .or. allocated(this%uniform_magnetic_field)
    
  end function hamiltonian_base_has_magnetic

  ! ----------------------------------------------------------------------------------

  logical pure function hamiltonian_base_projector_self_overlap(this) result(projector_self_overlap)
    type(hamiltonian_base_t), intent(in) :: this
    
    projector_self_overlap = this%projector_self_overlap
  end function hamiltonian_base_projector_self_overlap

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_base_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_base_inc.F90"

end module hamiltonian_base_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
