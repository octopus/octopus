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

module hamiltonian_elec_base_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use blas_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use global_oct_m
  use hardware_oct_m
  use hgh_projector_oct_m
  use ions_oct_m
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
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                                         &
    hamiltonian_elec_base_t,                        &
    dhamiltonian_elec_base_local,                   &
    zhamiltonian_elec_base_local,                   &
    dhamiltonian_elec_base_local_sub,               &
    zhamiltonian_elec_base_local_sub,               &
    dhamiltonian_elec_base_magnetic,                &
    zhamiltonian_elec_base_magnetic,                &
    dhamiltonian_elec_base_nlocal_start,            &
    zhamiltonian_elec_base_nlocal_start,            &
    dhamiltonian_elec_base_nlocal_finish,           &
    zhamiltonian_elec_base_nlocal_finish,           &
    dhamiltonian_elec_base_nlocal_position_commutator, &
    zhamiltonian_elec_base_nlocal_position_commutator, &
    hamiltonian_elec_base_has_magnetic,             &
    hamiltonian_elec_base_init,                     &
    hamiltonian_elec_base_end,                      &
    hamiltonian_elec_base_allocate,                 &
    hamiltonian_elec_base_clear,                    &
    hamiltonian_elec_base_build_proj,               &
    hamiltonian_elec_base_update,                   &
    hamiltonian_elec_base_accel_copy_pot,           &
    hamiltonian_elec_base_phase,                    &
    hamiltonian_elec_base_phase_spiral,             &
    hamiltonian_elec_base_rashba,                   &
    dhamiltonian_elec_base_nlocal_force,            &
    zhamiltonian_elec_base_nlocal_force,            &
    projection_t,                                   &
    hamiltonian_elec_base_projector_self_overlap,   &
    hamiltonian_elec_base_set_phase_corr,           &
    hamiltonian_elec_base_unset_phase_corr

  !> This object stores and applies an electromagnetic potential that
  !! can be represented by different types of potentials.

  type hamiltonian_elec_base_t
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
    logical,                               public :: has_non_local_potential
    integer                                       :: full_projection_size
    integer,                               public :: max_npoints
    integer,                               public :: total_points
    integer                                       :: max_nprojs
    logical                                       :: projector_mix
    CMPLX,                    allocatable, public :: projector_phases(:, :, :, :)
    integer,                  allocatable, public :: projector_to_atom(:)
    integer                                       :: nregions
    integer,                               public :: nphase
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
    CMPLX,                    allocatable, public :: phase(:, :)
    CMPLX,                    allocatable, public :: phase_corr(:,:)
    CMPLX,                    allocatable, public :: phase_spiral(:,:)
    type(accel_mem_t),                     public :: buff_phase
    type(accel_mem_t),                     public :: buff_phase_spiral
    integer,                               public :: buff_phase_qn_start
    logical                                       :: projector_self_overlap  !< if .true. some projectors overlap with themselves
    FLOAT,                    pointer,     public :: spin(:,:,:)
  end type hamiltonian_elec_base_t

  type projection_t
    private
    FLOAT, allocatable     :: dprojection(:, :)
    CMPLX, allocatable     :: zprojection(:, :)
    type(accel_mem_t)      :: buff_projection
    type(accel_mem_t)      :: buff_spin_to_phase
  end type projection_t

  integer, parameter, public ::          &
    TERM_ALL                 = HUGE(1),  &
    TERM_KINETIC             =   1,      &
    TERM_LOCAL_POTENTIAL     =   2,      & 
    TERM_NON_LOCAL_POTENTIAL =   4,      &
    TERM_OTHERS              =   8,      &
    TERM_LOCAL_EXTERNAL      =  16,      &
    TERM_MGGA                =  32,      &
    TERM_DFT_U               =  64,      &
    TERM_RDMFT_OCC           = 128

  integer, parameter, public ::            &
    FIELD_POTENTIAL                = 1,    &
    FIELD_VECTOR_POTENTIAL         = 2,    &
    FIELD_UNIFORM_VECTOR_POTENTIAL = 4,    &
    FIELD_UNIFORM_MAGNETIC_FIELD   = 8

  type(profile_t), save :: prof_vnlpsi_start, prof_vnlpsi_finish, prof_magnetic, prof_vlpsi, prof_scatter, &
    prof_matelement, prof_matelement_gather, prof_matelement_reduce

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_base_init(this, nspin, mass, rashba_coupling)
    type(hamiltonian_elec_base_t), intent(inout) :: this
    integer,                  intent(in)    :: nspin
    FLOAT,                    intent(in)    :: mass
    FLOAT,                    intent(in)    :: rashba_coupling

    PUSH_SUB(hamiltonian_elec_base_init)

    this%nspin = nspin
    this%mass  = mass
    this%rashba_coupling = rashba_coupling

    this%apply_projector_matrices = .false.
    this%has_non_local_potential = .false.
    this%nprojector_matrices = 0

    nullify(this%spin)
    this%projector_self_overlap = .false.
    
    POP_SUB(hamiltonian_elec_base_init)
  end subroutine hamiltonian_elec_base_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_base_end(this)
    type(hamiltonian_elec_base_t), intent(inout) :: this

    PUSH_SUB(hamiltonian_elec_base_end)

    if(allocated(this%potential) .and. accel_is_enabled()) then
      call accel_release_buffer(this%potential_opencl)
    end if
    
    SAFE_DEALLOCATE_A(this%potential)
    SAFE_DEALLOCATE_A(this%Impotential)
    SAFE_DEALLOCATE_A(this%vector_potential)
    SAFE_DEALLOCATE_A(this%uniform_vector_potential)
    SAFE_DEALLOCATE_A(this%uniform_magnetic_field)
    call hamiltonian_elec_base_destroy_proj(this)

    nullify(this%spin)

    POP_SUB(hamiltonian_elec_base_end)
  end subroutine hamiltonian_elec_base_end

  ! ---------------------------------------------------------- 
  !
  !> This functions sets to zero all fields that are currently
  !! allocated.
  !
  subroutine hamiltonian_elec_base_clear(this)
    type(hamiltonian_elec_base_t), intent(inout) :: this

    PUSH_SUB(hamiltonian_elec_clear)

    if(allocated(this%potential))                this%potential = M_ZERO
    if(allocated(this%Impotential))              this%Impotential = M_ZERO
    if(allocated(this%uniform_vector_potential)) this%uniform_vector_potential = M_ZERO
    if(allocated(this%vector_potential))         this%vector_potential = M_ZERO
    if(allocated(this%uniform_magnetic_field))   this%uniform_magnetic_field = M_ZERO

    POP_SUB(hamiltonian_elec_clear)
  end subroutine hamiltonian_elec_base_clear


  ! ---------------------------------------------------------------
  !> This function ensures that the corresponding field is allocated.
  subroutine hamiltonian_elec_base_allocate(this, mesh, field, complex_potential)
    type(hamiltonian_elec_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: field
    logical,                  intent(in)    :: complex_potential

    PUSH_SUB(hamiltonian_elec_base_allocate)

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

    POP_SUB(hamiltonian_elec_base_allocate)
  end subroutine hamiltonian_elec_base_allocate

  ! ---------------------------------------------------------- 
  !
  !> If both a uniform and non-uniform vector potentials are allocated,
  !! this function copies the uniform in the non-uniform one. In the
  !! future it may perform other internal consistency operations.
  !
  subroutine hamiltonian_elec_base_update(this, mesh)
    type(hamiltonian_elec_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh

    integer :: idir, ip

    PUSH_SUB(hamiltonian_elec_base_update)

    if(allocated(this%uniform_vector_potential) .and. allocated(this%vector_potential)) then
      ! copy the uniform vector potential onto the non-uniform one
      do idir = 1, mesh%sb%dim
        !$omp parallel do schedule(static)
        do ip = 1, mesh%np
          this%vector_potential(idir, ip) = &
            this%vector_potential(idir, ip) + this%uniform_vector_potential(idir)
        end do
      end do
      SAFE_DEALLOCATE_A(this%uniform_vector_potential)
    end if

    POP_SUB(hamiltonian_elec_base_update)
  end subroutine hamiltonian_elec_base_update


  !--------------------------------------------------------

  subroutine hamiltonian_elec_base_accel_copy_pot(this, mesh)
    type(hamiltonian_elec_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    
    integer :: offset, ispin

    PUSH_SUB(hamiltonian_elec_base_accel_copy_pot)

    if(allocated(this%potential) .and. accel_is_enabled()) then
      offset = 0
      do ispin = 1, this%nspin
        call accel_write_buffer(this%potential_opencl, mesh%np, this%potential(:, ispin), offset = offset)
        offset = offset + accel_padded_size(mesh%np)
      end do
    end if

    POP_SUB(hamiltonian_elec_base_accel_copy_pot)
  end subroutine hamiltonian_elec_base_accel_copy_pot

  
  !--------------------------------------------------------

  subroutine hamiltonian_elec_base_destroy_proj(this)
    type(hamiltonian_elec_base_t), intent(inout) :: this

    integer :: iproj

    PUSH_SUB(hamiltonian_elec_base_destroy_proj)

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

    POP_SUB(hamiltonian_elec_base_destroy_proj)
  end subroutine hamiltonian_elec_base_destroy_proj

  !-----------------------------------------------------------------
    
  subroutine hamiltonian_elec_base_build_proj(this, space, mesh, epot)
    type(hamiltonian_elec_base_t), target, intent(inout) :: this
    type(space_t),                         intent(in)    :: space
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

    PUSH_SUB(hamiltonian_elec_base_build_proj)

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
            overlap = submesh_overlap(epot%proj(iatom)%sphere, epot%proj(katom)%sphere, space)
            if(overlap) exit
          end do
        end if

        if(.not. overlap) then
          region_count(nregion) = region_count(nregion) + 1
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
          ASSERT(.not. submesh_overlap(epot%proj(order(iatom))%sphere, epot%proj(order(jatom))%sphere, space))
        end do
      end do
    end do

    call profiling_out(color_prof)

    ! deallocate previous projectors
    call hamiltonian_elec_base_destroy_proj(this)

    ! count projectors
    this%nprojector_matrices = 0
    this%apply_projector_matrices = .false.
    this%has_non_local_potential = .false.
    this%nregions = nregion

    !We determine if we have only local potential or not.
    do iorder = 1, epot%natoms
      iatom = order(iorder)

      if(.not. projector_is_null(epot%proj(iatom))) then
        this%has_non_local_potential = .true.
        exit
      end if
    end do

    do iorder = 1, epot%natoms
      iatom = order(iorder)
  
      if(projector_is(epot%proj(iatom), PROJ_KB) .or. projector_is(epot%proj(iatom), PROJ_HGH)) then
        this%nprojector_matrices = this%nprojector_matrices + 1
        this%apply_projector_matrices = .true.
        !The HGH pseudopotentials are now supporting the SOC
        if(epot%reltype /= NOREL .and. &
             (.not. projector_is(epot%proj(iatom), PROJ_HGH) .or. accel_is_enabled())) then 
          this%apply_projector_matrices = .false.
          exit
        end if
      else if(.not. projector_is_null(epot%proj(iatom))) then
        this%apply_projector_matrices = .false.
        exit
      end if
    end do

    if(mesh%use_curvilinear)  this%apply_projector_matrices = .false.

    if(.not. this%apply_projector_matrices) then
      SAFE_DEALLOCATE_A(order)
      SAFE_DEALLOCATE_A(head)

      POP_SUB(hamiltonian_elec_base_build_proj)
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
          
        iproj = iproj + 1

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
              nmat = nmat + epot%proj(iatom)%kb_p(ll, mm)%n_c
            end do
          end do
          
          call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%np, nmat, has_mix_matrix = .false.)
          
          ! generate the matrix
          pmat%dprojectors = M_ZERO
          imat = 1
          do ll = 0, lmax
            if (ll == lloc) cycle
            do mm = -ll, ll
              kb_p =>  epot%proj(iatom)%kb_p(ll, mm)
              do ic = 1, kb_p%n_c
                do ip = 1, pmat%npoints
                  pmat%dprojectors(ip, imat) = kb_p%p(ip, ic)
                end do
                pmat%scal(imat) = kb_p%e(ic)*mesh%vol_pp(1)
                imat = imat + 1
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
          
          call projector_matrix_allocate(pmat, epot%proj(iatom)%sphere%np, nmat, has_mix_matrix = .true., &
                                            is_cmplx = (epot%reltype == SPIN_ORBIT) )

          ! generate the matrix
          if(epot%reltype == SPIN_ORBIT) then
            pmat%zprojectors = M_ZERO
            pmat%zmix = M_ZERO
          else
            pmat%dprojectors = M_ZERO
            pmat%dmix = M_ZERO
          end if

          imat = 1
          do ll = 0, lmax
            if (ll == lloc) cycle
            do mm = -ll, ll
              hgh_p =>  epot%proj(iatom)%hgh_p(ll, mm)

              ! HGH pseudos mix different components, so we need to
              ! generate a matrix that mixes the projections
              if(epot%reltype == SPIN_ORBIT) then
                do ic = 1, 3
                  do jc = 1, 3
                    pmat%zmix(imat - 1 + ic, imat - 1 + jc, 1) = hgh_p%h(ic, jc) + M_HALF*mm*hgh_p%k(ic, jc)
                    pmat%zmix(imat - 1 + ic, imat - 1 + jc, 2) = hgh_p%h(ic, jc) - M_HALF*mm*hgh_p%k(ic, jc)
     
                    if(mm < ll) then
                      pmat%zmix(imat - 1 + ic, imat + 3 - 1 + jc, 3) = M_HALF*hgh_p%k(ic, jc)*sqrt(TOFLOAT(ll*(ll+1)-mm*(mm+1)))
                    end if

                    if(-mm < ll) then
                      pmat%zmix(imat - 1 + ic, imat - 3 - 1 + jc, 4) = M_HALF*hgh_p%k(ic, jc)*sqrt(TOFLOAT(ll*(ll+1)-mm*(mm-1)))
                    end if
                  end do
                end do
              else
                do ic = 1, 3
                  do jc = 1, 3
                    pmat%dmix(imat - 1 + ic, imat - 1 + jc) = hgh_p%h(ic, jc)
                  end do
                end do
              end if

              do ic = 1, 3
                if(epot%reltype == SPIN_ORBIT) then
                  do ip = 1, pmat%npoints
                    pmat%zprojectors(ip, imat) = hgh_p%zp(ip, ic)
                  end do
                else
                  do ip = 1, pmat%npoints
                    pmat%dprojectors(ip, imat) = hgh_p%dp(ip, ic)
                  end do
                end if
                pmat%scal(imat) = mesh%volume_element
                imat = imat + 1
              end do

            end do
          end do

          this%projector_self_overlap = this%projector_self_overlap .or. epot%proj(iatom)%sphere%overlap

        else
          cycle
        end if

        do ip = 1, pmat%npoints
          pmat%map(ip) = epot%proj(iatom)%sphere%map(ip)
          pmat%position(1:3, ip) = epot%proj(iatom)%sphere%x(ip, 1:3)
        end do

        this%full_projection_size = this%full_projection_size + pmat%nprojs

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
      this%total_points = this%total_points + pmat%npoints
    end do

    if(accel_is_enabled()) call build_opencl()

    POP_SUB(hamiltonian_elec_base_build_proj)

  contains

    subroutine build_opencl()
      integer              :: matrix_size, scal_size
      integer, allocatable :: cnt(:), invmap(:, :), invmap2(:), pos(:)
      integer, allocatable :: offsets(:, :)
      integer, parameter   :: OFFSET_SIZE = 6 ! also defined in share/opencl/projectors.cl
      integer, parameter   :: POINTS = 1, PROJS = 2, MATRIX = 3, MAP = 4, SCAL = 5, MIX = 6 ! update OFFSET_SIZE
      integer              :: ip, is, ii, ipos, mix_offset

      PUSH_SUB(hamiltonian_elec_base_build_proj.build_opencl)

      SAFE_ALLOCATE(offsets(1:OFFSET_SIZE, 1:this%nprojector_matrices))
      SAFE_ALLOCATE(cnt(1:mesh%np))

      cnt = 0

      ! Here we construct the offsets for accessing various arrays within the GPU kernels.
      ! The offset(:,:) array contains a number of sizes and offsets, describing how to address the arrays.
      ! This allows to transfer all these number to the GPU in one memory transfer.
      !
      ! For each projection matrix (addressed by imap), we have:
      !
      ! offset(POINTS, imap) : number of points of the sphere imap
      ! offset(PROJS, imap)  : number of projectors for imap
      ! offset(MATRIX, imap) : address offset: cumulative of pmat%npoints * pmat%nprojs
      ! offset(MAP, imap)    : address offset: cumulative of pmat%npoints for each imap
      ! offset(SCAL, imap)   : address_offset: cumulative of pmat%nprojs
      ! offset(MIX, imap)    : address_offset: cumulative of pmat%nprojs**2

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
        matrix_size = matrix_size + pmat%npoints*pmat%nprojs

        offsets(MAP, imat) = this%total_points
        this%total_points = this%total_points + pmat%npoints

        offsets(SCAL, imat) = scal_size
        scal_size = scal_size + pmat%nprojs

        if(allocated(pmat%dmix) .or. allocated(pmat%zmix)) then
          offsets(MIX, imat) = mix_offset
          mix_offset = mix_offset + pmat%nprojs**2
        else
          offsets(MIX, imat) = -1
        end if
        
        do is = 1, pmat%npoints
          ip = pmat%map(is)
          cnt(ip) = cnt(ip) + 1
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
          cnt(ip) = cnt(ip) + 1
          invmap(cnt(ip), ip) = ii
          ii = ii + 1
        end do
      end do

      ipos = 0
      pos(1) = 0
      do ip = 1, mesh%np
        do ii = 1, cnt(ip)
          ipos = ipos + 1
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
          call accel_write_buffer(this%buff_matrices, pmat%nprojs*pmat%npoints, pmat%dprojectors, offset = offsets(MATRIX, imat))
          call accel_write_buffer(this%buff_maps, pmat%npoints, pmat%map, offset = offsets(MAP, imat))
          call accel_write_buffer(this%buff_position, 3*pmat%npoints, pmat%position, offset = 3*offsets(MAP, imat))
        end if
        call accel_write_buffer(this%buff_scals, pmat%nprojs, pmat%scal, offset = offsets(SCAL, imat))
        if(offsets(MIX, imat) /= -1) call accel_write_buffer(this%buff_mix, pmat%nprojs**2, pmat%dmix, offset = offsets(MIX, imat))
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

      POP_SUB(hamiltonian_elec_base_build_proj.build_opencl)
    end subroutine build_opencl

  end subroutine hamiltonian_elec_base_build_proj
    
  ! ----------------------------------------------------------------------------------

  logical pure function hamiltonian_elec_base_has_magnetic(this) result(has_magnetic)
    type(hamiltonian_elec_base_t), intent(in) :: this
    
    has_magnetic = allocated(this%vector_potential) &
      .or. allocated(this%uniform_magnetic_field)
    
  end function hamiltonian_elec_base_has_magnetic

  ! ----------------------------------------------------------------------------------

  logical pure function hamiltonian_elec_base_projector_self_overlap(this) result(projector_self_overlap)
    type(hamiltonian_elec_base_t), intent(in) :: this
    
    projector_self_overlap = this%projector_self_overlap
  end function hamiltonian_elec_base_projector_self_overlap

 ! ----------------------------------------------------------------------------------

  subroutine hamiltonian_elec_base_set_phase_corr(hm_base, mesh, psib)
    type(hamiltonian_elec_base_t), intent(in) :: hm_base
    type(mesh_t),                  intent(in) :: mesh
    type(wfs_elec_t),           intent(inout) :: psib

    logical :: phase_correction

    PUSH_SUB(hamiltonian_elec_base_set_phase_corr)

    ! check if we only want a phase correction for the boundary points
    phase_correction = allocated(hm_base%phase)

    !We apply the phase only to np points, and the phase for the np+1 to np_part points
    !will be treated as a phase correction in the Hamiltonian
    if(phase_correction) then
      call hamiltonian_elec_base_phase(hm_base, mesh, mesh%np, .false., psib)
    end if

    POP_SUB(hamiltonian_elec_base_set_phase_corr)
  end subroutine hamiltonian_elec_base_set_phase_corr

  ! ----------------------------------------------------------------------------------

  subroutine hamiltonian_elec_base_unset_phase_corr(hm_base, mesh, psib)
    type(hamiltonian_elec_base_t), intent(in) :: hm_base
    type(mesh_t),                  intent(in) :: mesh
    type(wfs_elec_t),           intent(inout) :: psib

    logical :: phase_correction

    PUSH_SUB(hamiltonian_elec_base_unset_phase_corr)

    ! check if we only want a phase correction for the boundary points
    phase_correction = allocated(hm_base%phase)

    !We apply the phase only to np points, and the phase for the np+1 to np_part points
    !will be treated as a phase correction in the Hamiltonian
    if (phase_correction) then
      call hamiltonian_elec_base_phase(hm_base, mesh, mesh%np, .true., psib)
    end if

    POP_SUB(hamiltonian_elec_base_unset_phase_corr)
  end subroutine hamiltonian_elec_base_unset_phase_corr

  ! ---------------------------------------------------------------------------------------

  subroutine hamiltonian_elec_base_phase(this, mesh, np, conjugate, psib, src)
    type(hamiltonian_elec_base_t),         intent(in)    :: this
    type(mesh_t),                          intent(in)    :: mesh
    integer,                               intent(in)    :: np
    logical,                               intent(in)    :: conjugate
    type(wfs_elec_t),              target, intent(inout) :: psib
    type(wfs_elec_t),    optional, target, intent(in)    :: src

    integer :: ip, ii
    type(wfs_elec_t), pointer :: src_
    type(profile_t), save :: phase_prof
    CMPLX :: phase
    integer :: wgsize
    type(accel_kernel_t), save :: ker_phase

    PUSH_SUB(hamiltonian_elec_base_phase)
    call profiling_in(phase_prof, "PBC_PHASE_APPLY")

    call profiling_count_operations(6*np*psib%nst_linear)

    ASSERT(np <= mesh%np_part)
    ASSERT(psib%type() == TYPE_CMPLX)

    src_ => psib
    if(present(src)) src_ => src

    ASSERT(src_%has_phase .eqv. conjugate)
    ASSERT(src_%ik == psib%ik)
    ASSERT(src_%type() == TYPE_CMPLX)

    select case(psib%status())
    case(BATCH_PACKED)

      if(conjugate) then

        !$omp parallel do simd private(ii, phase)
        do ip = 1, np
          phase = conjg(this%phase(ip, psib%ik))
          do ii = 1, psib%nst_linear
            psib%zff_pack(ii, ip) = phase*src_%zff_pack(ii, ip)
          end do
        end do

      else

        !$omp parallel do simd private(ii, phase)
        do ip = 1, np
          phase = this%phase(ip, psib%ik)
          do ii = 1, psib%nst_linear
            psib%zff_pack(ii, ip) = phase*src_%zff_pack(ii, ip)
          end do
        end do

      end if

    case(BATCH_NOT_PACKED)

      if(conjugate) then

        !$omp parallel private(ii, ip)
        do ii = 1, psib%nst_linear
          !$omp do simd
          do ip = 1, np
            psib%zff_linear(ip, ii) = conjg(this%phase(ip, psib%ik))*src_%zff_linear(ip, ii)
          end do
          !$omp end do simd nowait
        end do
        !$omp end parallel

      else
        !$omp parallel private(ii, ip)
        do ii = 1, psib%nst_linear
          !$omp do simd
          do ip = 1, np
            psib%zff_linear(ip, ii) = this%phase(ip, psib%ik)*src_%zff_linear(ip, ii)
          end do
          !$omp end do simd nowait
        end do
        !$omp end parallel

      end if

    case(BATCH_DEVICE_PACKED)
      call accel_kernel_start_call(ker_phase, 'phase.cl', 'phase_hamiltonian')

      if(conjugate) then
        call accel_set_kernel_arg(ker_phase, 0, 1_4)
      else
        call accel_set_kernel_arg(ker_phase, 0, 0_4)
      end if

      call accel_set_kernel_arg(ker_phase, 1, (psib%ik - this%buff_phase_qn_start)*mesh%np_part)
      call accel_set_kernel_arg(ker_phase, 2, np)
      call accel_set_kernel_arg(ker_phase, 3, this%buff_phase)
      call accel_set_kernel_arg(ker_phase, 4, src_%ff_device)
      call accel_set_kernel_arg(ker_phase, 5, log2(src_%pack_size(1)))
      call accel_set_kernel_arg(ker_phase, 6, psib%ff_device)
      call accel_set_kernel_arg(ker_phase, 7, log2(psib%pack_size(1)))

      wgsize = accel_kernel_workgroup_size(ker_phase)/psib%pack_size(1)

      call accel_kernel_run(ker_phase, (/psib%pack_size(1), pad(np, wgsize)/), (/psib%pack_size(1), wgsize/))

      call accel_finish()
    end select

    psib%has_phase = .not. conjugate

    call profiling_out(phase_prof)
    POP_SUB(hamiltonian_elec_base_phase)
  end subroutine hamiltonian_elec_base_phase

  ! ---------------------------------------------------------------------------------------

  subroutine hamiltonian_elec_base_phase_spiral(this, der, psib)
    type(hamiltonian_elec_base_t),         intent(in)    :: this
    type(derivatives_t),                   intent(in)    :: der
    class(wfs_elec_t),                     intent(inout) :: psib

    integer               :: ip, ii, sp
    integer, allocatable  :: spin_label(:)
    type(accel_mem_t)     :: spin_label_buffer 
    type(profile_t), save :: phase_prof
    integer :: wgsize

    PUSH_SUB(hamiltonian_elec_base_phase_spiral)
    call profiling_in(phase_prof, "PBC_PHASE_SPIRAL")

    call profiling_count_operations(6*(der%mesh%np_part-der%mesh%np)*psib%nst_linear)

    ASSERT(der%boundaries%spiral)
    ASSERT(psib%type() == TYPE_CMPLX)

    sp = der%mesh%np
    if(der%mesh%parallel_in_domains) sp = der%mesh%np + der%mesh%vp%np_ghost


    select case(psib%status())
    case(BATCH_PACKED)

      !$omp parallel do private(ip, ii)
      do ip = sp + 1, der%mesh%np_part
        do ii = 1, psib%nst_linear, 2
          if(this%spin(3,psib%linear_to_ist(ii), psib%ik)>0) then
            psib%zff_pack(ii+1, ip) = psib%zff_pack(ii+1, ip)*this%phase_spiral(ip-sp, 1)
          else
            psib%zff_pack(ii, ip) = psib%zff_pack(ii, ip)*this%phase_spiral(ip-sp, 2)
          end if
        end do
      end do
      !$omp end parallel do

    case(BATCH_NOT_PACKED)

      !$omp parallel private(ii, ip)
      do ii = 1, psib%nst_linear, 2
        if(this%spin(3,psib%linear_to_ist(ii), psib%ik)>0) then
          !$omp do
          do ip = sp + 1, der%mesh%np_part
            psib%zff_linear(ip, ii+1) = psib%zff_linear(ip, ii+1)*this%phase_spiral(ip-sp, 1)
          end do
          !$omp end do nowait
        else
          !$omp do
          do ip = sp + 1, der%mesh%np_part
            psib%zff_linear(ip, ii) = psib%zff_linear(ip, ii)*this%phase_spiral(ip-sp, 2)
          end do
          !$omp end do nowait
        end if
      end do
      !$omp end parallel

    case(BATCH_DEVICE_PACKED)

      ASSERT(accel_is_enabled())

      ! generate array of offsets for access of psib and phase_spiral:
      ! TODO: Move this to the routine where spin(:,:,:) is generated
      !       and also move the buffer to the GPU at this point to
      !       avoid unecessary latency here!
 
      SAFE_ALLOCATE(spin_label(1:psib%nst_linear))
      spin_label = 0
      do ii = 1, psib%nst_linear, 2
        if(this%spin(3, psib%linear_to_ist(ii), psib%ik) > 0) spin_label(ii)=1
      end do

      call accel_create_buffer(spin_label_buffer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, psib%nst_linear)
      call accel_write_buffer(spin_label_buffer, psib%nst_linear, spin_label)

      call accel_kernel_start_call(kernel_phase_spiral, 'phase_spiral.cl', 'phase_spiral_apply')

      call accel_set_kernel_arg(kernel_phase_spiral, 0, psib%nst) 
      call accel_set_kernel_arg(kernel_phase_spiral, 1, sp) 
      call accel_set_kernel_arg(kernel_phase_spiral, 2, der%mesh%np_part) 
      call accel_set_kernel_arg(kernel_phase_spiral, 3, psib%ff_device)
      call accel_set_kernel_arg(kernel_phase_spiral, 4, log2(psib%pack_size(1)))
      call accel_set_kernel_arg(kernel_phase_spiral, 5, this%buff_phase_spiral)
      call accel_set_kernel_arg(kernel_phase_spiral, 6, spin_label_buffer)

      wgsize = accel_kernel_workgroup_size(kernel_phase_spiral)/psib%pack_size(1)

      call accel_kernel_run(kernel_phase_spiral, &
                            (/psib%pack_size(1)/2, pad(der%mesh%np_part - sp, 2*wgsize)/), &
                            (/psib%pack_size(1)/2, 2*wgsize/))

      call accel_finish()

      call accel_release_buffer(spin_label_buffer)

      SAFE_DEALLOCATE_A(spin_label)

    end select

    call profiling_out(phase_prof)
    POP_SUB(hamiltonian_elec_base_phase_spiral)
  end subroutine hamiltonian_elec_base_phase_spiral

  ! ---------------------------------------------------------------------------------------
  subroutine hamiltonian_elec_base_rashba(this, mesh, der, std, psib, vpsib)
    type(hamiltonian_elec_base_t),  intent(in)    :: this
    type(mesh_t),                   intent(in)    :: mesh
    type(derivatives_t),            intent(in)    :: der
    type(states_elec_dim_t),        intent(in)    :: std
    type(wfs_elec_t), target,       intent(in)    :: psib
    type(wfs_elec_t), target,       intent(inout) :: vpsib

    integer :: ist, idim, ip
    CMPLX, allocatable :: psi(:, :), vpsi(:, :), grad(:, :, :)

    PUSH_SUB(hamiltonian_elec_base_rashba)

    if(abs(this%rashba_coupling) < M_EPSILON) then
      POP_SUB(hamiltonian_elec_base_rashba)
      return
    end if
    ASSERT(std%ispin == SPINORS)
    ASSERT(mesh%sb%dim == 2)
    ASSERT(psib%type() == TYPE_CMPLX)
    ASSERT(vpsib%type() == TYPE_CMPLX)

    SAFE_ALLOCATE(psi(1:mesh%np_part, 1:std%dim))
    SAFE_ALLOCATE(vpsi(1:mesh%np, 1:std%dim))
    SAFE_ALLOCATE(grad(1:mesh%np, 1:mesh%sb%dim, 1:std%dim))

    do ist = 1, psib%nst
      call batch_get_state(psib, ist, mesh%np_part, psi)
      call batch_get_state(vpsib, ist, mesh%np, vpsi)

      do idim = 1, std%dim
        call zderivatives_grad(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
      end do

      if(allocated(this%vector_potential)) then
        do ip = 1, mesh%np
          vpsi(ip, 1) = vpsi(ip, 1) + &
            (this%rashba_coupling) * (this%vector_potential(2, ip) + M_zI * this%vector_potential(1, ip)) * psi(ip, 2)
          vpsi(ip, 2) = vpsi(ip, 2) + &
            (this%rashba_coupling) * (this%vector_potential(2, ip) - M_zI * this%vector_potential(1, ip)) * psi(ip, 1)
        end do
      end if

      do ip = 1, mesh%np
        vpsi(ip, 1) = vpsi(ip, 1) - &
          this%rashba_coupling*( grad(ip, 1, 2) - M_zI*grad(ip, 2, 2) )
        vpsi(ip, 2) = vpsi(ip, 2) + &
          this%rashba_coupling*( grad(ip, 1, 1) + M_zI*grad(ip, 2, 1) )
      end do

      call batch_set_state(vpsib, ist, mesh%np, vpsi)
    end do

    SAFE_DEALLOCATE_A(grad)
    SAFE_DEALLOCATE_A(vpsi)
    SAFE_DEALLOCATE_A(psi)
  
    POP_SUB(hamiltonian_elec_base_rashba)
  end subroutine hamiltonian_elec_base_rashba

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_elec_base_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_elec_base_inc.F90"

end module hamiltonian_elec_base_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
